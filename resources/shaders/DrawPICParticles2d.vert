#version 330 core
layout(location = 0) in vec2 pointPosition;

uniform float domainX;
uniform float domainY;
uniform float pointSize;

void main()
{
    // 将 [0, domain] 映射到 NDC [-1, 1]
    float x = (pointPosition.x / domainX) * 2.0 - 1.0;
    float y = (pointPosition.y / domainY) * 2.0 - 1.0;
    gl_Position = vec4(x, y, 0.0, 1.0);
    gl_PointSize = pointSize;
}
