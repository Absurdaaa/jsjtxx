#version 330 core

in vec2 vUV;
out vec4 FragColor;

uniform sampler3D densityTex;

uniform mat4 invViewProj;
uniform vec3 cameraPos;

uniform vec3 boxMin;
uniform vec3 boxMax;

uniform vec2 viewport;
uniform float stepSize;
uniform float densityScale;
uniform vec3 smokeColor;
uniform vec3 background;

bool intersectAABB(vec3 ro, vec3 rd, vec3 bmin, vec3 bmax, out float t0, out float t1)
{
    vec3 invD = 1.0 / rd;
    vec3 tBot = (bmin - ro) * invD;
    vec3 tTop = (bmax - ro) * invD;
    vec3 tMin = min(tBot, tTop);
    vec3 tMax = max(tBot, tTop);
    t0 = max(max(tMin.x, tMin.y), tMin.z);
    t1 = min(min(tMax.x, tMax.y), tMax.z);
    return t1 > max(t0, 0.0);
}

void main()
{
    // Reconstruct a world-space ray from this pixel.
    // OpenGL NDC: x,y in [-1,1], z=-1 near, z=+1 far.
    vec2 ndc = vUV * 2.0 - 1.0;

    vec4 nearH = invViewProj * vec4(ndc, -1.0, 1.0);
    vec4 farH  = invViewProj * vec4(ndc,  1.0, 1.0);
    vec3 nearP = nearH.xyz / nearH.w;
    vec3 farP  = farH.xyz / farH.w;

    vec3 ro = cameraPos;
    vec3 rd = normalize(farP - ro);

    float tEnter, tExit;
    if (!intersectAABB(ro, rd, boxMin, boxMax, tEnter, tExit))
    {
        FragColor = vec4(background, 1.0);
        return;
    }

    float t = max(tEnter, 0.0);
    vec3 boxSize = boxMax - boxMin;

    // Front-to-back compositing
    vec3 accumRGB = vec3(0.0);
    float accumA = 0.0;

    // Safety cap to avoid infinite loops
    int maxSteps = 1024;
    int steps = int((tExit - t) / stepSize) + 1;
    steps = clamp(steps, 1, maxSteps);

    for (int s = 0; s < steps; ++s)
    {
        vec3 p = ro + rd * t;
        vec3 uvw = (p - boxMin) / boxSize; // [0,1]

        float d = texture(densityTex, uvw).r;
        d = max(d, 0.0);

        // Convert density -> opacity per step (Beer-Lambert-ish)
        float a = 1.0 - exp(-d * densityScale * stepSize);
        a = clamp(a, 0.0, 1.0);

        vec3 c = smokeColor;
        // Optional: make very thin smoke dimmer
        c *= d;

        float oneMinusA = 1.0 - accumA;
        accumRGB += oneMinusA * a * c;
        accumA   += oneMinusA * a;

        if (accumA > 0.99)
            break;

        t += stepSize;
        if (t > tExit)
            break;
    }

    vec3 outRGB = accumRGB + (1.0 - accumA) * background;
    FragColor = vec4(outRGB, 1.0);
}
