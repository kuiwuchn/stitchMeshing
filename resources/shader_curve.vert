#version 330

in vec3 position;
in vec4 color;

uniform mat4 mvp;
out vec4 p;
out vec4 col_frag;

void main() {
	gl_Position = mvp * vec4(position, 1.0);
	col_frag = color;
	p = vec4(position, 1.0);
}
