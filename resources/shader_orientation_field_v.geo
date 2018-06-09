#version 330

layout(points) in;
layout(line_strip, max_vertices = 6) out;

in vData {
   vec3 q; 
   vec3 n;
} vertices[];

out vec4 color;

uniform mat4 mvp;
uniform float scale;

void main() {
    vec4 x = vec4(vertices[0].q, 0);
    vec4 y = vec4(cross(vertices[0].n, vertices[0].q), 0);

    vec4 p = gl_in[0].gl_Position;

	

color = vec4(1.0 , 1, 0.0, 1.0);
 gl_Position = mvp * (p - scale * x);
// gl_Position = mvp * (p);
   EmitVertex();
   gl_Position = mvp * (p + scale * x);
   EmitVertex();
   EndPrimitive();

color = vec4(0, 0, 1.0, 1.0);
   gl_Position = mvp * (p);
   EmitVertex();
   gl_Position = mvp * (p + scale *vec4(vertices[0].n, 0));
   EmitVertex();
   EndPrimitive();

//color = vec4(0 , 1, 1, 1.0);
//   gl_Position = mvp * (p - scale *y);
//   EmitVertex();
//   gl_Position = mvp * (p + scale *y);
//   EmitVertex();
//   EndPrimitive();
   
   
//    color = vec4(1.0, 0.0, 0.0, 1.0);
//    gl_Position = mvp * (p - scale*x);
//    EmitVertex();
//    gl_Position = mvp * (p + scale*x);
//    EmitVertex();
//    EndPrimitive();
//
//    color = vec4(0.0, 1.0, 0.0, 1.0);
//    gl_Position = mvp * (p - scale*y);
//    EmitVertex();
//    gl_Position = mvp * (p + scale*y);
//    EmitVertex();
//    EndPrimitive();
}
