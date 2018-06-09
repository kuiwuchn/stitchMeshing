#version 330

uniform mat4 mvp;

in vec3 o;
out vec4 p;

void main() {
    p = vec4(o, 1.0);
	gl_Position = mvp * vec4(o, 1.0);
}
