#version 330

in vec3 position;
in vec3 q;
in vec3 n;

out vData {
   vec3 q;
   vec3 n;
} geo;

void main() {
	gl_Position = vec4(position, 1.0);
	geo.q = q;
	geo.n = n;
}
