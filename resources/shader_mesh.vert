#version 330

in vec3 position;
in vec3 normal;

out vec3 normal_geo;

out vData {
    vec3 normal;
} geo;

void main() {
	gl_Position = vec4(position, 1.0);
	geo.normal = normal;
}
