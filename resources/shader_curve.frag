#version 330
precision lowp float;

out vec4 outColor;

uniform vec4 split;
in vec4 p;
in vec4 col_frag;

void main() {
    if (dot(split, p) < 0)
        discard;

    outColor = col_frag;
}
