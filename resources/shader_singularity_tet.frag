#version 330
precision lowp float;

out vec4 outColor;

uniform vec4 color;
uniform vec4 split;

in fData {
    vec3 color;
    vec4 p;
} frag;

void main() {
    if (dot(split, frag.p) < 0)
        discard;

    outColor = vec4(frag.color, 1.0);
}
