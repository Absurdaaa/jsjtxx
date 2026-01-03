#version 330 core

in float vDensity;
out vec4 FragColor;

void main()
{
    // 点精灵局部坐标 [0,1] → [-0.5,0.5]
    vec2 uv = gl_PointCoord - vec2(0.5);
    float r2 = dot(uv, uv);

    // 高斯核（烟雾关键）
    float alpha = exp(-r2 * 16.0);

    // 密度调制透明度
    alpha *= clamp(vDensity, 0.0, 1.0);

    // 边缘裁剪，防止方块
    if (alpha < 0.01)
        discard;

    // 烟雾颜色（浅灰）
    vec3 smokeColor = vec3(0.8);

    FragColor = vec4(smokeColor, alpha);
}
