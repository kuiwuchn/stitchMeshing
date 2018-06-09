#version 330

in vec3 position;

uniform mat4 mvp;
out vec4 p;

void main() {
	gl_Position = mvp * vec4(position, 1.0);
	p = vec4(position, 1.0);
}
