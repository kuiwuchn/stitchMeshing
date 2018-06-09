#version 330
precision lowp float;

out vec4 outColor;

uniform vec4 split;
in vec4 p;

void main() {
    if (dot(split, p) < 0)
        discard;

    outColor = vec4(0.0, 0.0, 1.0, 0.1);
}
