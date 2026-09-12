#ifndef FUZZY_PID_H
#define FUZZY_PID_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ===================================================================
 * 模糊自调度 PID（只调度 Kp / Kd，Ki 固定 + 反算抗饱和）
 *
 * 坐标约定（务必保持）：
 *   err = 目标位置 - 小球当前位置 = 0 - pos
 *   pos 向左为负、向右为正，管道范围 [-10, +10]
 *   ec  = d(err)/dt = -小球速度
 * 因此：
 *   球在左   -> err > 0
 *   球在右   -> err < 0
 *   球向右滚 -> ec  < 0
 *   球向左滚 -> ec  > 0
 * =================================================================== */

/* ---------- 物理量满量程（决定量化因子 KE / KEC） ---------- */
#ifndef FUZZY_E_MAX
#define FUZZY_E_MAX       10.0f
#endif

#ifndef FUZZY_EC_MAX
#define FUZZY_EC_MAX      5.0f
#endif

/* ---------- PID 基准增益 ---------- */
#ifndef FUZZY_KP_BASE
#define FUZZY_KP_BASE     1.0f
#endif
#ifndef FUZZY_KD_BASE
#define FUZZY_KD_BASE     0.5f
#endif
#ifndef FUZZY_KI_FIXED
#define FUZZY_KI_FIXED    0.4f
#endif

/* ---------- 模糊增量 -> 实际增益的换算因子 ---------- */
#ifndef FUZZY_OUT_KP
#define FUZZY_OUT_KP      0.05f
#endif
#ifndef FUZZY_OUT_KD
#define FUZZY_OUT_KD      0.05f
#endif

/* ---------- 增益限幅 ---------- */
#ifndef FUZZY_KP_MIN
#define FUZZY_KP_MIN      0.5f
#endif
#ifndef FUZZY_KP_MAX
#define FUZZY_KP_MAX      3.0f
#endif
#ifndef FUZZY_KD_MIN
#define FUZZY_KD_MIN      0.2f
#endif
#ifndef FUZZY_KD_MAX
#define FUZZY_KD_MAX      2.0f
#endif

/* ---------- 输出限幅（控制量 u，具体单位由执行器换算） ---------- */
#ifndef FUZZY_OUT_MIN
#define FUZZY_OUT_MIN    -30.0f
#endif
#ifndef FUZZY_OUT_MAX
#define FUZZY_OUT_MAX     30.0f
#endif

/* ---------- 反算抗饱和系数（越大，退出饱和越快） ---------- */
#ifndef FUZZY_KAW
#define FUZZY_KAW         0.2f
#endif

/* ---------- 积分安全兜底 ----------
 * 注意：这不是正常抗饱和手段，正常抗饱和由 back-calculation 完成。
 * 这里只设一个很宽的边界，防止数值发散/溢出；正常工作时不应碰到它。
 */
#ifndef FUZZY_I_SAFETY
#define FUZZY_I_SAFETY    (2.0f * FUZZY_OUT_MAX)
#endif

/* ---------- dt 保护 ---------- */
#ifndef FUZZY_DT_MIN
#define FUZZY_DT_MIN      1e-4f
#endif
#ifndef FUZZY_DT_MAX
#define FUZZY_DT_MAX      0.1f
#endif

/* ---------- EC 一阶低通滤波系数 ---------- */
#ifndef FUZZY_LP_EC
#define FUZZY_LP_EC       0.2f
#endif

/* ---------- 输入误差保护 ---------- */
#ifndef FUZZY_ERR_LIMIT
#define FUZZY_ERR_LIMIT   20.0f
#endif

typedef struct {
    float prev_err;    /* 上一次误差 */
    float ec_f;        /* 低通滤波后的误差变化率 */
    float integral;    /* 积分器状态 */
    float kp;
    float kd;
    float ki;
    float out;         /* 上一次输出，异常时保持 */
    uint8_t first_call;
} FuzzyPID;

void fuzzy_pid_init(FuzzyPID *pid);
float fuzzy_pid_step(FuzzyPID *pid, float err, float dt);

#ifdef __cplusplus
}
#endif

#endif /* FUZZY_PID_H */