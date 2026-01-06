#version 330 core
out vec4 FragColor;

in vec2 TexCoord;

uniform sampler2D densityTex;  // 密度场纹理
uniform float contrast;        // 对比度

// --- 小工具：hash噪声(用于轻微抖动，降低条带感) ---
float hash12(vec2 p)
{
    // 低成本 hash，足够做dither
    vec3 p3  = fract(vec3(p.xyx) * 0.1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}

// Catmull-Rom cubic
float cubicCR(float v0, float v1, float v2, float v3, float t)
{
    float a0 = -0.5*v0 + 1.5*v1 - 1.5*v2 + 0.5*v3;
    float a1 = v0 - 2.5*v1 + 2.0*v2 - 0.5*v3;
    float a2 = -0.5*v0 + 0.5*v2;
    float a3 = v1;
    return ((a0*t + a1)*t + a2)*t + a3;
}

// 16-tap bicubic（在低分辨率密度图上比GL_LINEAR更细腻）
float sampleBicubic(sampler2D tex, vec2 uv)
{
    ivec2 texSizeI = textureSize(tex, 0);
    vec2 texSize = vec2(texSizeI);
    vec2 invTex = 1.0 / texSize;

    // 像素空间坐标（以像素中心对齐）
    vec2 xy = uv * texSize - 0.5;
    vec2 ixy = floor(xy);
    vec2 f = xy - ixy;

    // 取 4x4（GLSL 330 不允许“数组的数组”，用一维数组展开）
    float s[16];
    for (int j = 0; j < 4; ++j)
    {
        for (int i = 0; i < 4; ++i)
        {
            vec2 o = (ixy + vec2(i - 1, j - 1) + 0.5) * invTex;
            s[j * 4 + i] = texture(tex, o).r;
        }
    }

    float col0 = cubicCR(s[0 * 4 + 0], s[0 * 4 + 1], s[0 * 4 + 2], s[0 * 4 + 3], f.x);
    float col1 = cubicCR(s[1 * 4 + 0], s[1 * 4 + 1], s[1 * 4 + 2], s[1 * 4 + 3], f.x);
    float col2 = cubicCR(s[2 * 4 + 0], s[2 * 4 + 1], s[2 * 4 + 2], s[2 * 4 + 3], f.x);
    float col3 = cubicCR(s[3 * 4 + 0], s[3 * 4 + 1], s[3 * 4 + 2], s[3 * 4 + 3], f.x);
    return cubicCR(col0, col1, col2, col3, f.y);
}

void main()
{
    // 采样密度（bicubic重建）
    float d = sampleBicubic(densityTex, TexCoord);
    d = max(d, 0.0);

    // 对比度（保持你现有UI参数）
    d = pow(d, 1.0 / max(contrast, 0.1));

    // 指数映射：把“密度动态范围”压到 [0,1)，同时保持层次
    // k 越大越“厚/更不透明”（推荐 0.2~1.0）
    const float k = 0.55;
    float dVis = 1.0 - exp(-k * d);

    // 指数吸收：alpha 用 dVis 驱动会更稳定（避免 d 稍大就几乎全白）
    const float sigma = 1.0;
    float alpha = 1.0 - exp(-sigma * dVis);

    // 轻微抖动：减少banding（幅度很小，基本看不见噪点）
    float n = hash12(gl_FragCoord.xy);
    alpha = clamp(alpha + (n - 0.5) * 0.01, 0.0, 1.0);

    // 基于密度梯度的“伪法线”做简单光照
    ivec2 ts = textureSize(densityTex, 0);
    vec2 texel = 1.0 / vec2(ts);
    float dx = sampleBicubic(densityTex, TexCoord + vec2(texel.x, 0.0)) - sampleBicubic(densityTex, TexCoord - vec2(texel.x, 0.0));
    float dy = sampleBicubic(densityTex, TexCoord + vec2(0.0, texel.y)) - sampleBicubic(densityTex, TexCoord - vec2(0.0, texel.y));
    vec2 grad = vec2(dx, dy);
    vec2 N = normalize(vec2(-grad.x, -grad.y) + 1e-6);
    vec2 L = normalize(vec2(-0.55, 0.85)); // 光方向：从左下照向右上（可按喜好调整）
    float lambert = clamp(0.5 + 0.5 * dot(N, L), 0.0, 1.0);

    // 简单自阴影：沿光反方向采样几步累积厚度
    // 步数/强度越大越“立体”，但也更耗。
    float shadow = 1.0;
    {
        const int steps = 10;
        const float stepLen = 2.0; // 每步跨几个 texel
        float acc = 0.0;
        vec2 stepUV = (-L) * texel * stepLen;
        vec2 uv = TexCoord;
        for (int i = 0; i < steps; ++i)
        {
            uv += stepUV;
            float sd = sampleBicubic(densityTex, uv);
            acc += max(sd, 0.0);
        }
        shadow = exp(-0.18 * acc);
    }

    // 颜色：冷色烟（可按喜好换色/做color ramp）
    vec3 bgColor = vec3(0.02, 0.02, 0.03);
    vec3 smokeDark = vec3(0.10, 0.11, 0.13);
    vec3 smokeLight = vec3(0.88, 0.90, 0.94);

    vec3 smoke = mix(smokeDark, smokeLight, clamp(dVis, 0.0, 1.0));
    smoke *= mix(0.75, 1.15, lambert) * shadow;

    vec3 color = mix(bgColor, smoke, alpha);
    FragColor = vec4(color, 1.0);
}
