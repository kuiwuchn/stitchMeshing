#version 330

in vec3 position;
in vec3 color;

uniform mat4 mvp;

out fData {
    vec3 color;
    vec4 p;
} frag;

void main() {
    gl_Position = mvp * vec4(position, 1.0);
    frag.p = vec4(position, 1.0);
    frag.color = color;
}
