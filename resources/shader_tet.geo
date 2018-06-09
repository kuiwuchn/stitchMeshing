#version 330

layout(lines_adjacency) in;
layout(triangle_strip, max_vertices = 12) out;

uniform vec3 light_position;
uniform mat4 proj, model, view;
uniform vec4 split;

in vData {
    vec4 color;
} vertices[];

out fData {
    vec3 to_eye;
    vec3 to_light;
    vec3 normal;
    vec4 color;
} frag;

void main() {
    int indices[12] = int[](
       1, 0, 2,
       3, 2, 0,
       1, 2, 3,
       0, 1, 3
    );

    if (dot(split, gl_in[0].gl_Position) < 0 ||
        dot(split, gl_in[1].gl_Position) < 0 ||
        dot(split, gl_in[2].gl_Position) < 0 ||
        dot(split, gl_in[3].gl_Position) < 0)
        return;

    for (int i=0; i<4; ++i) {
        int id[3] = int[](
           indices[3*i+0],
           indices[3*i+1],
           indices[3*i+2]
        );
        vec3 normal = normalize(cross(
            gl_in[id[1]].gl_Position.xyz-gl_in[id[0]].gl_Position.xyz,
            gl_in[id[2]].gl_Position.xyz-gl_in[id[0]].gl_Position.xyz));
        frag.normal = (model * (view * vec4(normal, 0.0))).xyz;

        for (int j=0; j<3; ++j) {
            vec4 pos_camera = view * (model * gl_in[id[j]].gl_Position);
            gl_Position = proj * pos_camera;
            frag.to_light = (view * vec4(light_position, 1.0)).xyz - pos_camera.xyz;
            frag.to_eye = -pos_camera.xyz;
            frag.color = vertices[id[j]].color;
            EmitVertex();
        }

        EndPrimitive();
    }
}
