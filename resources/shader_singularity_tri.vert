#version 330

in vec3 position;
in vec3 normal;
in vec3 color;

out vData {
    vec3 color;
    vec3 normal;
} vertex;

void main() {
    gl_Position = vec4(position, 1.0);
    vertex.normal = normal;
    vertex.color = color;
}
