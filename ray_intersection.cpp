#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <random>
using namespace std;

inline double random_double() {
    return rand() / (RAND_MAX + 1.0);
}

struct vec3{
    double x, y, z;

    vec3() : x(0), y(0), z(0) {}
    vec3(double e0, double e1, double e2) : x(e0), y(e1), z(e2) {}

    double length() const {
        return sqrt(x*x + y*y + z*z);
    }
    
    double length_squared() const {
        return x*x + y*y + z*z;
    }
};

struct material{
    vec3 color;   // RGB 
    double ambient; // [0, 1]
    double specular; //[0, 1]
    double phong; //Power
    double refractive_index;

    material() : color(0, 0, 0), ambient(0.0), specular(0.0), phong(1), refractive_index(1) {}
    material(vec3 a, double b, double c, double d, double e) : color(a.x, a.y, a.z), ambient(b), specular(c), phong(d), refractive_index(e) {}
};

struct hit_record{
    double t;
    vec3 p, normal;
    material mat;

    hit_record() : t(0), p(0, 0, 0), normal(0, 0, 0), mat(vec3(0, 0, 0), 0, 0, 1, 1.0) {}
    hit_record(double a1, vec3 a2, vec3 a3, material m) : t(a1), p(a2.x, a2.y, a2.z), normal(a3.x, a3.y, a3.z), mat(m.color, m.ambient, m.specular, m.phong, m.refractive_index) {}
};


inline vec3 operator+(const vec3 &u, const vec3 &v) {
    return vec3(u.x + v.x, u.y + v.y, u.z + v.z);
}

inline vec3 operator-(const vec3 &u, const vec3 &v) {
    return vec3(u.x - v.x, u.y - v.y, u.z - v.z);
}

inline vec3 operator*( const double d, const vec3 &u) {
    return vec3(u.x * d , u.y * d, u.z * d);
}

inline vec3 operator*(const vec3 &u, const double d) {
    return vec3(u.x * d , u.y * d, u.z * d);
}

inline vec3 operator/(const vec3 &u, const double d) {
    return vec3(u.x / d , u.y / d, u.z / d);
}

inline vec3 normalize(const vec3 &u) {
    double val = sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
    return u / val;
}

inline double vec_length(const vec3 &u) {
    double val = sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
    return val;
}

inline vec3 operator*(const vec3 &u, const vec3 &v) {
    return vec3(u.x * v.x, u.y * v.y, u.z * v.z);
}

inline double dot(const vec3 &u, const vec3 &v) {
    return u.x * v.x
         + u.y * v.y
         + u.z * v.z;
}

inline vec3 cross(const vec3 &u, const vec3 &v) {
    return vec3(u.y * v.z - u.z * v.y,
                u.z * v.x - u.x * v.z,
                u.x * v.y - u.y * v.x);
}

//For soft shadows -> Nudged light pos should be inside sphere
inline vec3 random_in_unit_sphere() {
    while (true) {

        // random_double() * 2.0 - 1.0 gives a range of [-1.0 to 1.0]
        vec3 p = vec3(random_double()*2.0 - 1.0, 
                      random_double()*2.0 - 1.0, 
                      random_double()*2.0 - 1.0);
        
        // If it's inside the sphere -> keep it Otherwise, the loop runs again.
        if (p.length_squared() < 1.0) return p;
    }
}

vec3 refracted_ray(const vec3& r_incomming, const hit_record& rec) {
    vec3 uv = normalize(r_incomming); 
    double dt = dot(uv, rec.normal);
    
    double discriminant;
    double eta; //n_i / n_t
    vec3 outward_normal;

    if (dt > 0) { //Inside to Outside
        outward_normal = -1 * rec.normal;  
        eta = rec.mat.refractive_index / 1.0; // Material / Air
        dt = dot(uv, outward_normal); 
    } 

    else { //Outside to Inside
        outward_normal = rec.normal;    // Normal is already correct
        eta = 1.0 / rec.mat.refractive_index; // Air / Glass
    }

    // Discriminant in Snell's Law :1.0 - eta^2 * (1 - cos^2(theta))
    discriminant = 1.0 - eta * eta * (1.0 - dt * dt);

    if (discriminant < 0) { //Total Internal Reflectionn
        // Return a standard reflection vector
        return uv - 2 * dot(uv, outward_normal) * outward_normal;
    }

    vec3 refracted = eta * (uv - outward_normal * dt) - outward_normal * sqrt(discriminant);
    return refracted;
}

double schlick(double cosine, double ref_idx) {
    double r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1.0 - r0) * pow((1.0 - cosine), 5);
}

class Surface {
    public:
    virtual bool hit(const vec3& r_origin, const vec3& r_dir, double t0, double t1, hit_record& rec) = 0;
};

class Sphere: public Surface{
    public:
        vec3 center;
        double radius;
        material mat;
    
    Sphere(vec3 c, double r, material m) : center(c), radius(r), mat(m) {}

    //Checking if a ray from eye(e + td) will pass through a sphere (|x - c| = r)
    virtual bool hit(const vec3& r_origin, const vec3& r_dir, double t0, double t1, hit_record& rec) override {
        vec3 oc = r_origin - center; 
        double a = dot(r_dir, r_dir);
        double b = 2.0 * dot(oc, r_dir);
        double c_val = dot(oc, oc) - radius*radius; 
        double discriminant = b*b - 4*a*c_val;

        if (discriminant < 0) return false;

        double sqrtd = sqrt(discriminant);
        double root = (-b - sqrtd) / (2.0 * a);
        
        if (root < t0 || root > t1) {
            root = (-b + sqrtd) / (2.0 * a);
            if (root < t0 || root > t1) return false;
        }

        rec.t = root;
        rec.p = r_origin + rec.t * r_dir;
        rec.normal = (rec.p - center) / radius;
        rec.mat = this->mat;
        return true;
    }
};

class Triangle: public Surface{
    public:
        vec3 a, b, c;
        material mat;

    Triangle(vec3 v0, vec3 v1, vec3 v2, material m) : a(v0), b(v1), c(v2), mat(m) {}

    // Checking if a ray from eye(e) with direction(d) will intersect triangle (a,b,c)
    virtual bool hit(const vec3& r_origin, const vec3& r_dir, double t0, double t1, hit_record& rec) override {
        double A = a.x - b.x;
        double B = a.y - b.y;
        double C = a.z - b.z; 
        
        double D = a.x - c.x;
        double E = a.y - c.y;
        double F = a.z - c.z; 
        
        double G = r_dir.x;
        double H = r_dir.y;
        double I = r_dir.z;
        
        double J = a.x - r_origin.x;
        double K = a.y - r_origin.y;
        double L = a.z - r_origin.z;

        double EI_HF = E * I - H * F;
        double GF_DI = G * F - D * I;
        double DH_EG = D * H - E * G;

        double M = A * (EI_HF) + B * (GF_DI) + C * (DH_EG);

        if (abs(M) < 1e-8) return false; 

        double AK_JB = A * K - J * B;
        double JC_AL = J * C - A * L;
        double BL_KC = B * L - K * C;

        double t = - (F * AK_JB + E * JC_AL + D * BL_KC) / M;
        
        if (t < t0 || t > t1) {
            return false;
        }

        double gam = (I * AK_JB + H * JC_AL + G * BL_KC) / M;
        
        if (gam < 0 || gam > 1) {
            return false;
        }

        double bet = (J * EI_HF + K * GF_DI + L * DH_EG) / M;

        if (bet < 0 || bet > 1 - gam) {
            return false;
        }

        rec.t = t;
        rec.p = r_origin + rec.t * r_dir;
        rec.normal = normalize(cross(b-a, c-a));
        rec.mat = this->mat;

        return true;
    }
};

class Plane: public Surface{
    public:
        vec3 centre, normal;
        material mat;
        bool background;

    Plane(vec3 a, vec3 b, material m, bool c): centre(a), normal(b), mat(m), background(c) {}

    virtual bool hit(const vec3& r_origin, const vec3& r_dir, double t0, double t1, hit_record& rec) override {
        double den = dot(r_dir, normal);
        if(abs(den) < 1e-8) return false;
        
        double t = dot(centre - r_origin, normal)/den;
        if (t < t0 || t > t1) {
            return false;
        }

        rec.t = t;
        rec.p = r_origin + rec.t * r_dir;
        rec.normal = normal;
        rec.mat = this->mat;

        if(background){
            if((abs(int(floor(rec.p.x)) + int(floor(rec.p.z))) & 1) == 1){
                rec.mat.color = vec3(0.3, 0.3, 0.3);
            }
            else {
                rec.mat.color = vec3(0.9, 0.9, 0.9);
            }   
        }
        return true;
    }
};

class SurfaceList: public Surface{
    public:
        Surface** list; // Pointer to an array of Surface pointers
        int list_size;

    SurfaceList(Surface** l, int n) {
        list = l;
        list_size = n;
    }

    virtual bool hit(const vec3& r_origin, const vec3& r_dir, double t0, double t1, hit_record& rec) override {
            hit_record temp_rec;
            bool hit_anything = false;
            double closest_so_far = t1; // Start with the max distance (infinity)

            for (int i = 0; i < list_size; i++) {
                if (list[i]->hit(r_origin, r_dir, t0, closest_so_far, temp_rec)) {
                    hit_anything = true;
                    closest_so_far = temp_rec.t;
                    rec = temp_rec;
                }
            }
            return hit_anything;
        }
};

vec3 ray_tracing(const vec3& ray_origin, const vec3& ray_dir, const vec3& light_pos,  Surface* world, int depth){
    if (depth <= 0) {
        return vec3(0, 0, 0);
    }

    hit_record rec;

    if (world->hit(ray_origin, ray_dir, 0.001, 1000.0, rec)) {

        vec3 light_color(1, 1, 1); 
        
        vec3 l = normalize(light_pos - rec.p); // Direction of light
        vec3 v = normalize(ray_origin - rec.p); // Direction of eye
        vec3 h = normalize(l + v); // Halfway vector

        // Refraction
        if (rec.mat.refractive_index > 1.0) {
            vec3 refraction_color(0,0,0);
            vec3 reflection_color(0,0,0);
            
            vec3 refracted_dir;
            vec3 next_path = refracted_ray(ray_dir, rec);
            vec3 next_origin = rec.p + (next_path * 0.001);
            refraction_color = rec.mat.color * ray_tracing(next_origin, next_path, light_pos, world, depth - 1); //Multiplying with color for tint
            
            vec3 reflected_dir = ray_dir - (rec.normal * 2.0 * dot(ray_dir, rec.normal));
            vec3 reflected_origin = rec.p + (rec.normal * 0.001);
            reflection_color = ray_tracing(reflected_origin, reflected_dir, light_pos ,world, depth - 1);

            vec3 unit_dir = normalize(ray_dir);
            double cosine = dot(-1 * unit_dir, rec.normal);

            if (dot(unit_dir, rec.normal) > 0) {
                cosine = dot(unit_dir, rec.normal) * rec.mat.refractive_index; 
            }

            double R = schlick(cosine, rec.mat.refractive_index);

            vec3 final_glass = (refraction_color * (1.0 - R)) + (reflection_color * R);

            // Add the specular highlight (light bulb shine) on top
            double spec = pow(max(0.0, dot(rec.normal, h)), rec.mat.phong);
            vec3 specular_light = rec.mat.specular * spec * light_color;

            return final_glass + specular_light;
        }

        vec3 ambient = rec.mat.ambient * rec.mat.color;  //c = c_r * c_a

        vec3 shadow_ray_origin = rec.p + (rec.normal * 0.001); // To avoid self intersection
        double dist_to_light = (light_pos - rec.p).length(); //As light is not at infinity

        hit_record shadow_rec;
        vec3 color = ambient;

        vec3 shadow_attenuation(1.0, 1.0, 1.0);

        for( int i = 0; i < 10; i++)
        {
            bool in_shadow = world->hit(shadow_ray_origin, l, 0.001, dist_to_light, shadow_rec);
            //The light sorce should not cast shadow from reflection to it
            if(!in_shadow || shadow_rec.mat.ambient >= 1.0) break;

            //Opaque object -> Complete Black Shadow
            if(shadow_rec.mat.refractive_index <= 1.0){
                shadow_attenuation = vec3(0, 0, 0);
                break;
            }

            //For glasss -> Shadow will not be cast and instead light will go through with tint
            shadow_attenuation = shadow_attenuation * shadow_rec.mat.color;
            shadow_ray_origin = shadow_rec.p + l * 0.001;
            dist_to_light = (light_pos - shadow_ray_origin).length(); 
        }

        double diff = max(0.0, dot(rec.normal, l));
        vec3 diffuse = diff * rec.mat.color * light_color;

        double spec = pow(max(0.0, dot(rec.normal, h)), rec.mat.phong);
        vec3 specular = rec.mat.specular * spec * light_color;
        color = color + (diffuse * shadow_attenuation) + (specular * shadow_attenuation);

        if (rec.mat.specular > 0 && rec.mat.refractive_index <= 1.0) {
             vec3 reflected_dir = ray_dir - (rec.normal * 2.0 * dot(ray_dir, rec.normal));
             reflected_dir = normalize(reflected_dir);
             vec3 reflected_origin = rec.p + (rec.normal * 0.001);
             
             vec3 reflected_color = ray_tracing(reflected_origin, reflected_dir, light_pos, world, depth - 1);
             color = color + (reflected_color * rec.mat.specular);
        }

        return color;
    }
    //If ray does not hit -> Surface behind something
    return vec3(0.2, 0.2, 0.2);
}


//PATH TRACING
inline void create_coordinate_system(const vec3& N, vec3& Nt, vec3& Nb) {
    if (abs(N.x) > abs(N.y)) { 
        vec3 temp(N.z, 0.0, -N.x);
        Nt = normalize(temp);
    }
    else { 
        vec3 temp(0.0, -N.z, N.y);
        Nt = normalize(temp);
    }
    Nb = cross(N, Nt); 
}

vec3 uniform_sample_hemisphere(const float&r1, const float&r2) {
    float cos_theta = r1;
    float sin_theta = sqrt(1 - cos_theta * cos_theta);
    float phi = 2 * M_PI * r2;

    float x = sin_theta * cosf(phi);
    float z = sin_theta * sinf(phi);

    return vec3(x, r1, z); //y = cos_theta
}


vec3 path_tracing(const vec3& ray_origin, const vec3& ray_dir, Sphere* light_sphere,  Surface* world, int depth, bool first_bounce) {
    if (depth <= 0) {
        return vec3(0, 0, 0);
    }

    hit_record rec;

    if (world->hit(ray_origin, ray_dir, 0.001, 1000.0, rec)) {
        //Hitting a light source
        if (rec.mat.ambient >= 1.0) {
            if(first_bounce)
                return rec.mat.color; // Color of the light source
            else {
                return vec3(0, 0, 0);
            }
        }

        //Checking for Glass
        if (rec.mat.refractive_index > 1.0) {
            vec3 unit_dir = normalize(ray_dir);
            double cosine = dot(-1.0 * unit_dir, rec.normal);
            if (dot(unit_dir, rec.normal) > 0) cosine = dot(unit_dir, rec.normal) * rec.mat.refractive_index;
        
            double R = schlick(cosine, rec.mat.refractive_index);

            if (random_double() < R) {
            vec3 reflected_dir = normalize(ray_dir - (rec.normal * 2.0 * dot(ray_dir, rec.normal)));
            vec3 reflected_origin = rec.p + (rec.normal * 0.001);
            return path_tracing(reflected_origin, reflected_dir, light_sphere,  world, depth - 1, true);
            }
            else {
            vec3 refracted_dir = refracted_ray(ray_dir, rec);
            vec3 refracted_origin = rec.p + (refracted_dir * 0.001);
            return rec.mat.color * path_tracing(refracted_origin, refracted_dir, light_sphere, world, depth - 1, true); // Tinted refraction
            }
        }

        //Checking for Mirror/metal (Specular Reflection)
        if (rec.mat.specular > 0.0) {
            if (random_double() < rec.mat.specular) {
                vec3 reflected_dir = normalize(ray_dir - (rec.normal * 2.0 * dot(ray_dir, rec.normal)));
                vec3 reflected_origin = rec.p + (rec.normal * 0.001);
                return path_tracing(reflected_origin, reflected_dir, light_sphere, world, depth - 1, true); 
            }
        }

        //Base Case: Diffuse Reflection
        vec3 direct_light_color(0, 0, 0);
        vec3 indirect_light_color(0, 0, 0);


        vec3 light_target = light_sphere->center + (random_in_unit_sphere() * light_sphere->radius);
        vec3 light_dir = light_target - rec.p; 
        
        float dist_squared = light_dir.length_squared();
        float dist_to_light = sqrt(dist_squared);
        light_dir = normalize(light_dir);

        float r1 = random_double(); 
        float r2 = random_double();
        vec3 sample = uniform_sample_hemisphere(r1, r2);

        vec3 Nt, Nb;
        create_coordinate_system(rec.normal, Nt, Nb);

        vec3 sampleWorld( 
            sample.x * Nb.x + sample.y * rec.normal.x + sample.z * Nt.x,
            sample.x * Nb.y + sample.y * rec.normal.y + sample.z * Nt.y,
            sample.x * Nb.z + sample.y * rec.normal.z + sample.z * Nt.z
        );

        vec3 scattered_dir = normalize(sampleWorld);
        vec3 next_origin = rec.p + (rec.normal * 0.001); 

        // Direct Lighting
        hit_record shadow_rec;
        bool in_shadow = world->hit(next_origin, light_dir, 0.001, dist_to_light, shadow_rec);

        // If we hit the light source
        if(in_shadow && shadow_rec.mat.ambient >= 1.0){

            float cos_theta = max(0.0, dot(rec.normal, light_dir));
            float attenuation = 1.0 / dist_squared; 
            
            direct_light_color = rec.mat.color * light_sphere->mat.color * cos_theta * attenuation * 15.0;
        }

        //Indirect Lighting
        vec3 incoming_indirect_color = path_tracing(next_origin, scattered_dir, light_sphere, world, depth - 1, false);
        indirect_light_color = rec.mat.color * incoming_indirect_color * r1 * 2.0;

        return direct_light_color + indirect_light_color;
    } 
    return vec3(0.2, 0.2, 0.2);
}


int main() {
    // 1. Image Settings
    const int image_width = 400;
    const int image_height = 400;

    // 2. Camera Settings
    // The "Viewport" is a virtual screen floating in front of the camera
    double viewport_height = 2.0;
    double viewport_width = 2.0;
    double focal_length = 1.0;

    vec3 origin(0, 0, 0);
    vec3 horizontal(viewport_width, 0, 0);
    vec3 vertical(0, viewport_height, 0);
    
    // Calculate the position of the lower-left corner of the viewport
    // Origin - Half Width - Half Height - Focal Length (Depth)
    vec3 lower_left_corner = origin - (horizontal/2.0) - (vertical/2.0) - vec3(0, 0, focal_length);

    // 3. Define Materials
    material glass(vec3(1.0, 0.5, 0.5), 0.0, 1.0, 100.0, 1.5); //Some specular reflection with big value of phong 
    material red_plastic(vec3(0.7, 0.1, 0.1), 0.1, 0.0, 50.0, 1.0); 
    material green_matte(vec3(0.1, 0.8, 0.1), 0.5, 0.0, 1.0, 1.0);
    material blue_metal(vec3(0.2, 0.2, 0.8), 0.2, 0.8, 20.0, 1.0);
    material tile_grey(vec3(0.3, 0.3, 0.3), 0.2, 0.0, 20.0, 1.0);
    material light_source(vec3(5.0, 5.0, 5.0), 1.0, 0.0, 0.0, 1.0); 

    Surface* my_objects[5];

    // Object 0: The Target (Red Sphere)
    my_objects[0] = new Sphere(vec3(0, 0, -5), 1.0, glass); 
    
    // Object 1: The Background (Green Triangle)
    my_objects[1] = new Triangle(vec3(-3,-2,-6), vec3(1,-2,-6), vec3(-1,2,-6), green_matte);

    // Object 2: The Shadow Caster (Blue Sphere)
    // Placed between the light (10,10,10) and the Red Sphere to cast shadow and get reflections
    my_objects[2] = new Sphere(vec3(1.0, 0.5, -2.5), 0.5, blue_metal);

    //Objet 3: A sphere acting as the light source 
    //Base Light position -> Centre of Sphere
    vec3 base_light_pos(4.0, 4.0, -3.0);
    double base_light_radius = 2;
    Sphere light_sphere(base_light_pos, base_light_radius, light_source);
    my_objects[3] = new Sphere(base_light_pos, base_light_radius, light_source);

    //Object 4: The background tiles
    my_objects[4] = new Plane(vec3(0, -3, 0), vec3(0, 1, 0), tile_grey, true);

    SurfaceList world(my_objects, 5);
    
    const int samples_per_pixel = 200; //For Anti-Aliasing and Monte Carlo Estimation

    bool path = false;

    // 5. Render Loop (Outputting PPM format)
    // Header: P3 means ASCII color, then Width Height, then Max Color Value (255)
    cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    // Loop through every row (j) from Top to Bottom
    for (int j = image_height - 1; j >= 0; --j) {
        
        // Loop through every column (i) from Left to Right
        for (int i = 0; i < image_width; ++i) {

            vec3 pixel_color(0, 0, 0);
            vec3 indirect_diffuse(0, 0, 0);

            for (int k = 0; k < samples_per_pixel; k++) {

                // A. Calculate UV Coordinates (0.0 to 1.0)
                // This maps pixel (200, 200) to the center of the viewport (0.5, 0.5)
                double u = double(i + random_double()) / (image_width - 1);
                double v = double(j + random_double()) / (image_height - 1);

                vec3 sample_light_pos = base_light_pos + (random_in_unit_sphere() * base_light_radius);
                // B. Create the Ray
                // Start at the eye (origin) and go through the specific spot on the viewport
                vec3 ray_direction = lower_left_corner + (u * horizontal) + (v * vertical) - origin;
                
                // C. Calculate Color (The Physics Engine)
                if(!path)
                    pixel_color = pixel_color + ray_tracing(origin, normalize(ray_direction), sample_light_pos, &world, 5);
                else{
                    indirect_diffuse = indirect_diffuse + path_tracing(origin, normalize(ray_direction), &light_sphere, &world, 5, true);
                }
            }
            vec3 final_color;
            if(!path) {
                final_color = pixel_color / samples_per_pixel;
            } else {
                final_color = indirect_diffuse / samples_per_pixel; 
            }

            //Clamping
            final_color.x = final_color.x / (final_color.x + 1.0);
            final_color.y = final_color.y / (final_color.y + 1.0);
            final_color.z = final_color.z / (final_color.z + 1.0);

            // D. Write Color to Output
            //Gamma correction
            final_color.x = sqrt(final_color.x);
            final_color.y = sqrt(final_color.y);
            final_color.z = sqrt(final_color.z);
            
            // Convert 0.0-1.0 range to 0-255 integers
            cout << static_cast<int>(255.999 * final_color.x) << ' '
                 << static_cast<int>(255.999 * final_color.y) << ' '
                 << static_cast<int>(255.999 * final_color.z) << '\n';
        }
    }

    return 0;
}