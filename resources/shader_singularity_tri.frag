#version 330
precision lowp float;

in fData {
    vec3 color;
    vec2 texcoord;
} frag;

out vec4 outColor;

void main() {
    if (length(frag.texcoord.xy-vec2(0.5)) > 0.5)
        discard;
    vec3 col = frag.color;
    outColor = vec4(col, 1.0);
}
