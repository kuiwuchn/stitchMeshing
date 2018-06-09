#version 330

layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;

uniform vec3 light_position;
uniform mat4 proj, model, view;
uniform vec4 split;

in vData {
    vec3 normal;
} vertices[];

out fData {
	vec3 to_eye;
	vec3 to_light;
	vec3 normal;
	vec4 color;
} frag;

void main() {
if (dot(split, gl_in[0].gl_Position) < 0 ||
        dot(split, gl_in[1].gl_Position) < 0 ||
        dot(split, gl_in[2].gl_Position) < 0)
        return;
		
    frag.normal = normalize(cross(
        (view * (model * (gl_in[1].gl_Position - gl_in[0].gl_Position))).xyz,
        (view * (model * (gl_in[2].gl_Position - gl_in[0].gl_Position))).xyz));

    for (int i=0; i<3; ++i) {
        vec4 pos_camera = view * (model * gl_in[i].gl_Position);
        gl_Position = proj * pos_camera;
        frag.to_light = (view * vec4(light_position, 1.0)).xyz - pos_camera.xyz;
        frag.to_eye = -pos_camera.xyz;
        frag.color = vec4(0.0);//vertices[i].normal != vec3(0.0) ? vec4(0.8, 0.0, 0.0, 0.2) : vec4(0.0, 0.0, 1.0, 0.2);
        EmitVertex();
    }

    EndPrimitive();
}
