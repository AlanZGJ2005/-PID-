#include "fuzzy_pid.h"

/* ===================================================================
 * 模糊论域：[-3, +3]，7 档：NB NM NS ZO PS PM PB
 *
 * 量化因子：
 *   KE  = 3 / FUZZY_E_MAX
 *   KEC = 3 / FUZZY_EC_MAX
 *
 * 定点说明：
 *   论域坐标用 Q13（1 个模糊档间距 = 8192）
 *   隶属度/权重用 Q15（1.0 = 32768）
 * =================================================================== */

#define N_MF            7
#define Q13_ONE         8192
#define Q13_SCALE       8192.0f
#define Q13_LIMIT       24576

#define KE   (3.0f / FUZZY_E_MAX)
#define KEC  (3.0f / FUZZY_EC_MAX)

/* 7 个档：NB NM NS ZO PS PM PB */
enum {
    FUZZY_NB = 0,
    FUZZY_NM,
    FUZZY_NS,
    FUZZY_ZO,
    FUZZY_PS,
    FUZZY_PM,
    FUZZY_PB
};

/* 各档中心：-3, -2, -1, 0, 1, 2, 3（Q13） */
static const int32_t MF_CENTER_Q13[N_MF] = {
    -24576, -16384, -8192, 0, 8192, 16384, 24576
};

/* ===================================================================
 * 规则表（行 = E，列 = EC）
 *
 * 列/行顺序：NB NM NS ZO PS PM PB
 *
 * 表义：
 *   KP：球越远越加“拉力”；但球已经在快速往回冲时，要收一点力。
 *   KD：速度越快，阻尼越大；中心附近额外加一点，防止冲过头。
 *
 * 示例（按本文件 err/ec 约定）：
 *   球在 -8cm、正以 +6cm/s 向右回中：
 *     err = +8 -> E 在 PM/PB 之间
 *     ec  = -6 -> EC 落在 NB
 *     重点看 [PB][NB]：KP=1（收力）、KD=3（大阻尼）
 * =================================================================== */
static const int8_t KP_RULE[N_MF][N_MF] = {
    /*        EC: NB NM NS ZO PS PM PB */
    /* E NB */ { 3, 3, 3, 3, 2, 2, 1 },
    /* E NM */ { 3, 3, 3, 2, 1, 1, 0 },
    /* E NS */ { 3, 2, 2, 1, 0, 0,-1 },
    /* E ZO */ { 1, 1, 0, 0, 0, 1, 1 },
    /* E PS */ {-1, 0, 0, 1, 2, 2, 3 },
    /* E PM */ { 0, 1, 1, 2, 3, 3, 3 },
    /* E PB */ { 1, 2, 2, 3, 3, 3, 3 },
};

static const int8_t KD_RULE[N_MF][N_MF] = {
    /*        EC: NB NM NS ZO PS PM PB */
    /* E NB */ { 3, 2, 1,-1, 1, 2, 3 },
    /* E NM */ { 3, 2, 1,-1, 1, 2, 3 },
    /* E NS */ { 3, 3, 2, 0, 2, 3, 3 },
    /* E ZO */ { 3, 3, 2, 0, 2, 3, 3 },
    /* E PS */ { 3, 3, 2, 0, 2, 3, 3 },
    /* E PM */ { 3, 2, 1,-1, 1, 2, 3 },
    /* E PB */ { 3, 2, 1,-1, 1, 2, 3 },
};

/* ---------- 基础工具 ---------- */
static float clampf(float v, float lo, float hi)
{
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

/* 物理量 -> Q13 论域，并限幅、四舍五入 */
static int32_t scale_to_q13(float x, float k)
{
    float v = x * k * Q13_SCALE;

    if (v >= (float)Q13_LIMIT) {
        return Q13_LIMIT;
    }
    if (v <= -(float)Q13_LIMIT) {
        return -Q13_LIMIT;
    }

    if (v >= 0.0f) {
        return (int32_t)(v + 0.5f);
    }
    return (int32_t)(v - 0.5f);
}

/*
 * 找到 x 落在哪两个相邻三角隶属函数之间：
 *   idx      : 左侧档的下标
 *   mu_left  : 左侧档隶属度（Q15，1.0 = 32768）
 *   mu_right : 右侧档隶属度（Q15）
 */
static void mf_two_q15(int32_t x_q13,
                       int32_t *idx,
                       int32_t *mu_left,
                       int32_t *mu_right)
{
    /* x_q13 + Q13_LIMIT 一定非负，右移是确定行为 */
    int32_t i = (x_q13 + Q13_LIMIT) >> 13;
    if (i < 0) i = 0;
    if (i > N_MF - 1) i = N_MF - 1;

    int32_t frac = x_q13 - MF_CENTER_Q13[i];
    if (frac < 0) frac = 0;
    if (frac > Q13_ONE) frac = Q13_ONE;

    *idx = i;
    *mu_left  = (Q13_ONE - frac) << 2;  /* Q13 -> Q15 */
    *mu_right = frac << 2;
}

/* 累加一个规则格子的贡献 */
static void add_cell(int32_t ei,
                     int32_t eci,
                     int32_t mu_e,
                     int32_t mu_ec,
                     int32_t *num_kp,
                     int32_t *num_kd,
                     int32_t *den)
{
    int32_t w = (mu_e * mu_ec) >> 15;   /* Q15 权重 */

    if (w <= 0) {
        return;
    }

    *num_kp += w * (int32_t)KP_RULE[ei][eci];
    *num_kd += w * (int32_t)KD_RULE[ei][eci];
    *den    += w;
}

/*
 * 模糊推理：E、EC -> dkp、dkd
 * dkp/dkd 是 [-3, +3] 的模糊增量，还没有乘输出因子。
 */
static void fuzzy_gain_q15(float e, float ec, float *dkp, float *dkd)
{
    int32_t e_q13  = scale_to_q13(e, KE);
    int32_t ec_q13 = scale_to_q13(ec, KEC);

    int32_t ei, eci;
    int32_t meL, meR, mcL, mcR;
    int32_t num_kp = 0;
    int32_t num_kd = 0;
    int32_t den = 0;

    mf_two_q15(e_q13,  &ei,  &meL, &meR);
    mf_two_q15(ec_q13, &eci, &mcL, &mcR);

    /* 最多 2 x 2 = 4 个格子有贡献 */
    if (meL > 0) {
        if (mcL > 0) {
            add_cell(ei, eci, meL, mcL, &num_kp, &num_kd, &den);
        }
        if (mcR > 0 && eci < N_MF - 1) {
            add_cell(ei, eci + 1, meL, mcR, &num_kp, &num_kd, &den);
        }
    }

    if (meR > 0 && ei < N_MF - 1) {
        if (mcL > 0) {
            add_cell(ei + 1, eci, meR, mcL, &num_kp, &num_kd, &den);
        }
        if (mcR > 0 && eci < N_MF - 1) {
            add_cell(ei + 1, eci + 1, meR, mcR, &num_kp, &num_kd, &den);
        }
    }

    if (den > 0) {
        *dkp = (float)num_kp / (float)den;
        *dkd = (float)num_kd / (float)den;
    } else {
        *dkp = 0.0f;
        *dkd = 0.0f;
    }
}

/* ---------- 初始化 ---------- */
void fuzzy_pid_init(FuzzyPID *pid)
{
    if (pid == NULL) {
        return;
    }

    pid->prev_err = 0.0f;
    pid->ec_f = 0.0f;
    pid->integral = 0.0f;
    pid->kp = FUZZY_KP_BASE;
    pid->kd = FUZZY_KD_BASE;
    pid->ki = FUZZY_KI_FIXED;
    pid->out = 0.0f;
    pid->first_call = 1;
}

/* ---------- 单步计算 ---------- */
float fuzzy_pid_step(FuzzyPID *pid, float err, float dt)
{
    float ec_raw, dkp, dkd;
    float p, d, i;
    float u_raw, u;

    if (pid == NULL) {
        return 0.0f;
    }

    /* 1. 输入保护：NaN / 传感器异常 / dt 异常时，保持上一次输出 */
    if (err != err) {
        return pid->out;
    }
    if (err > FUZZY_ERR_LIMIT || err < -FUZZY_ERR_LIMIT) {
        return pid->out;
    }
    if (dt < FUZZY_DT_MIN) {
        return pid->out;
    }
    if (dt > FUZZY_DT_MAX) {
        dt = FUZZY_DT_MAX;
    }

    /* 2. 计算误差变化率；首次调用用 0，避免假大 EC 冲击 D */
    if (pid->first_call) {
        pid->prev_err = err;
        ec_raw = 0.0f;
        pid->first_call = 0;
    } else {
        ec_raw = (err - pid->prev_err) / dt;
        pid->prev_err = err;
    }

    /* 3. EC 一阶低通滤波，抑制 K230 像素噪声 */
    pid->ec_f += FUZZY_LP_EC * (ec_raw - pid->ec_f);

    /* 4. 模糊查表：得到 dkp、dkd */
    fuzzy_gain_q15(err, pid->ec_f, &dkp, &dkd);

    pid->kp = clampf(FUZZY_KP_BASE + FUZZY_OUT_KP * dkp, FUZZY_KP_MIN, FUZZY_KP_MAX);
    pid->kd = clampf(FUZZY_KD_BASE + FUZZY_OUT_KD * dkd, FUZZY_KD_MIN, FUZZY_KD_MAX);

    /* 5. PID 三项 */
    p = pid->kp * err;
    d = pid->kd * pid->ec_f;
    i = pid->integral;

    u_raw = p + i + d;
    u = clampf(u_raw, FUZZY_OUT_MIN, FUZZY_OUT_MAX);

    /* 6. 反算抗饱和：正常抗饱和主力 */
    pid->integral += pid->ki * err * dt;
    pid->integral += FUZZY_KAW * (u - u_raw) * dt;

    /* 7. 积分安全兜底：只防溢出/发散，正常工作时不应碰到 */
    pid->integral = clampf(pid->integral, -FUZZY_I_SAFETY, FUZZY_I_SAFETY);

    pid->out = u;
    return u;
}
