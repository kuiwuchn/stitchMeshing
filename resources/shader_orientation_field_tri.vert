#version 330

in vec3 position;
in vec3 q;
in vec3 n;
in vec3 scale_color;

out vData {
   vec3 q;
   vec3 n;
   vec3 scale_color;   
} geo;

void main() {
	gl_Position = vec4(position, 1.0);
	geo.q = q;
	geo.n = n;
geo.scale_color = scale_color;	
}
