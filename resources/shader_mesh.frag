#version 330
precision lowp float;

uniform vec4 specular_color;
uniform vec4 base_color;

in fData {
	vec3 to_eye;
	vec3 to_light;
	vec3 normal;
	vec4 color;
} frag;

out vec4 outColor;

void main() {
    vec4 Kd = base_color;
    vec4 Ks = specular_color;

    if (frag.color.a != 0.0) {
        Kd = frag.color;
    }
    vec4 Ka = Kd * 0.2;

    vec3 to_light = normalize(frag.to_light);
    vec3 to_eye = normalize(frag.to_eye);
    vec3 normal = normalize(frag.normal);
    vec3 refl = reflect(-to_light, normal);
    if (dot(to_light, normal) <= 0)
        discard;

    float diffuse_factor = max(0.0, dot(to_light, normal));
    float specular_factor = pow(max(dot(to_eye, refl), 0.0), 10.0);

    outColor = Ka + Kd*diffuse_factor + Ks*specular_factor;
}
