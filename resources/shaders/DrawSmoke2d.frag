#version 330 core
out vec4 FragColor;

in vec2 TexCoord;

uniform sampler2D densityTex;  // 密度场纹理
uniform float contrast;        // 对比度

void main()
{
    // 采样密度值
    float density = texture(densityTex, TexCoord).r;

    // 应用对比度
    density = pow(density, 1.0 / max(contrast, 0.1));

    // 烟雾颜色：从深灰到浅白
    vec3 smokeColor = mix(vec3(0.1, 0.1, 0.12), vec3(0.9, 0.92, 0.95), density);

    // 背景色（深色）
    vec3 bgColor = vec3(0.02, 0.02, 0.03);

    // 混合烟雾和背景
    vec3 finalColor = mix(bgColor, smokeColor, density);

    FragColor = vec4(finalColor, 1.0);
}
