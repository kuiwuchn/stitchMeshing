#version 330

layout(points) in;
layout(triangle_strip, max_vertices = 4) out;

uniform mat4 mvp;
uniform float scale;

in vData {
    vec3 color;
    vec3 normal;
} vertices[];

out fData {
    vec3 color;
    vec2 texcoord;
} frag;


const vec2 corners[4] = vec2[](vec2(0.0, 1.0), vec2(0.0, 0.0), vec2(1.0, 1.0), vec2(1.0, 0.0));

void main() {
    vec3 n = vertices[0].normal;
    vec3 s = vec3(1.0, 2.0, 4.5);
    s = normalize(s - n*dot(n, s));
    vec3 t = cross(n, s);

    for (int i=0; i<4; ++i) {
        vec4 pos = gl_in[0].gl_Position;
        vec2 corner = corners[i];
        vec3 d = s * (corner.x - 0.5) + t * (corner.y - 0.5);
        pos.xyz += scale * d;
        gl_Position = mvp * pos;
        frag.texcoord = corners[i];
        frag.color = vertices[0].color;
        EmitVertex();
    }
    EndPrimitive();
}
