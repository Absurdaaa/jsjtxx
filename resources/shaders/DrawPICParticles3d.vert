#version 330 core

// 粒子世界坐标
layout (location = 0) in vec3 aPos;

// （可选）粒子密度 / 权重
layout (location = 1) in float aDensity;

uniform mat4 view;
uniform mat4 projection;

// 基础点大小（像素）
uniform float basePointSize;

out float vDensity;

void main()
{
    vDensity = aDensity;

    // 世界 → 裁剪空间
    vec4 viewPos = view * vec4(aPos, 1.0);
    gl_Position = projection * viewPos;

    // 透视缩放（近大远小）
    float dist = length(viewPos.xyz);
    gl_PointSize = basePointSize / max(dist, 0.1);
}
