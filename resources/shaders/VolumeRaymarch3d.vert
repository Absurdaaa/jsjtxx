#version 330 core

// fullscreen triangle (no VBO). Outputs uv in [0,1].
out vec2 vUV;

void main()
{
    // 3 vertices: (-1,-1), (3,-1), (-1,3)
    vec2 pos;
    if (gl_VertexID == 0) pos = vec2(-1.0, -1.0);
    else if (gl_VertexID == 1) pos = vec2( 3.0, -1.0);
    else pos = vec2(-1.0,  3.0);

    vUV = pos * 0.5 + 0.5;
    gl_Position = vec4(pos, 0.0, 1.0);
}
