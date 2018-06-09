#version 330

in vec3 position;
in vec4 color;

out vData {
    vec4 color;
} geo;

void main() {
	gl_Position = vec4(position, 1.0);
	geo.color = color;
}
