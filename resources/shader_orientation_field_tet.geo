#version 330

layout(points) in;
layout(line_strip, max_vertices = 6) out;

in vData {
   vec4 q; 
   vec3 scale_color;
} vertices[];

out vec4 color;

uniform mat4 mvp;
uniform vec4 split;
uniform vec4 split2;
uniform float scale;

void main() {
    vec4 q = vertices[0].q;

    float xx = q.x * q.x, yy = q.y * q.y, zz = q.z * q.z;
    float xy = q.x * q.y, xz = q.x * q.z, yz = q.y * q.z;
    float xw = q.x * q.w, yw = q.y * q.w, zw = q.z * q.w;

    vec4 x = vec4(1 - 2*yy - 2*zz, 2*(xy + zw), 2*(xz - yw), 0.0);
    vec4 y = vec4(2*(xy - zw), 1 - 2*xx - 2*zz, 2*(yz + xw), 0.0);
    vec4 z = vec4(2*(xz + yw), 2*(yz - xw), 1 - 2*xx - 2*yy, 0.0);

    vec4 p = gl_in[0].gl_Position;

    if (dot(split, gl_in[0].gl_Position) < 0)
        return;
    if (dot(split2, gl_in[0].gl_Position) > 0)
    return;	

//    color = vec4(1.0, 0.0, 0.0, 1.0);
//    gl_Position = mvp * (p - vertices[0].scale[0]*x);
//    EmitVertex();
//    gl_Position = mvp * (p + vertices[0].scale[0]*x);
//    EmitVertex();
//    EndPrimitive();
//
//    color = vec4(0.0, 1.0, 0.0, 1.0);
//    gl_Position = mvp * (p - vertices[0].scale[1]*y);
//    EmitVertex();
//    gl_Position = mvp * (p + vertices[0].scale[1]*y);
//    EmitVertex();
//    EndPrimitive();
//
//    color = vec4(0.0, 0.0, 1.0, 1.0);
//    gl_Position = mvp * (p - vertices[0].scale[2]*z);
//    EmitVertex();
//    gl_Position = mvp * (p + vertices[0].scale[2]*z);
//    EmitVertex();
//    EndPrimitive();

	if(vertices[0].scale_color[0]>=-0.0000000001 &&vertices[0].scale_color[0]<=0.0000000001)
	color = vec4(1.0 * vertices[0].scale_color[0], 1- vertices[0].scale_color[0], 0.0, 0.0);
else color = vec4(1.0 * vertices[0].scale_color[0], 1- vertices[0].scale_color[0], 0.0, 1.0);
   gl_Position = mvp * (p - scale * x);
   EmitVertex();
   gl_Position = mvp * (p + scale * x);
   EmitVertex();
   EndPrimitive();

if(vertices[0].scale_color[1]>=-0.0000000001 &&vertices[0].scale_color[1]<=0.0000000001)
	color = vec4(1.0 * vertices[0].scale_color[1], 1- vertices[0].scale_color[1], 0.0, 0.0);
else color = vec4(1.0 * vertices[0].scale_color[1], 1- vertices[0].scale_color[1], 0.0, 1.0);
   gl_Position = mvp * (p - scale *y);
   EmitVertex();
   gl_Position = mvp * (p + scale *y);
   EmitVertex();
   EndPrimitive();

if(vertices[0].scale_color[2]>=-0.0000000001 &&vertices[0].scale_color[2]<=0.0000000001)
	color = vec4(1.0 * vertices[0].scale_color[2], 1- vertices[0].scale_color[2], 0.0, 0.0);
else color = vec4(1.0 * vertices[0].scale_color[2], 1- vertices[0].scale_color[2], 0.0, 1.0);
   gl_Position = mvp * (p - scale *z);
   EmitVertex();
   gl_Position = mvp * (p + scale *z);
   EmitVertex();
   EndPrimitive();

//color = vec4(1.0, 0.0, 0.0, 1.0);
//gl_Position = mvp * (p - scale*x);
//EmitVertex();
//gl_Position = mvp * (p + scale*x);
//EmitVertex();
//EndPrimitive();
//
//color = vec4(0.0, 1.0, 0.0, 1.0);
//gl_Position = mvp * (p - scale*y);
//EmitVertex();
//gl_Position = mvp * (p + scale*y);
//EmitVertex();
//EndPrimitive();
//
//color = vec4(0.0, 0.0, 1.0, 1.0);
//gl_Position = mvp * (p - scale*z);
//EmitVertex();
//gl_Position = mvp * (p + scale*z);
//EmitVertex();
//EndPrimitive();
}
