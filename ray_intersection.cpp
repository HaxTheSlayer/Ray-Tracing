#include <iostream>
#include <cmath>
using namespace std;

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
    double phong;   

    material() : color(0, 0, 0), ambient(0.0), specular(0.0), phong(1) {}
    material(vec3 a, double b, double c, double d) : color(a.x, a.y, a.z), ambient(b), specular(c), phong(d) {}
};

struct hit_record{
    double t;
    vec3 p, normal;
    material mat;

    hit_record() : t(0), p(0, 0, 0), normal(0, 0, 0), mat(vec3(0, 0, 0), 0, 0, 1) {}
    hit_record(double a1, vec3 a2, vec3 a3, material m) : t(a1), p(a2.x, a2.y, a2.z), normal(a3.x, a3.y, a3.z), mat(m.color, m.ambient, m.specular, m.phong) {}
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

vec3 hit_color(const vec3& ray_origin, const vec3& ray_dir, Surface* world){
    hit_record rec;

    if (world->hit(ray_origin, ray_dir, 0.001, 1000.0, rec)) {

        vec3 light_pos(10, 10, 10); 
        vec3 light_color(1, 1, 1); 
        
        vec3 l = normalize(light_pos - rec.p); // Direction of light
        vec3 v = normalize(ray_origin - rec.p); // Direction of eye
        vec3 h = normalize(l + v); // Halfway vector

        vec3 ambient = rec.mat.ambient * rec.mat.color;  //c = c_r * c_a

        vec3 shadow_ray_origin = rec.p + (rec.normal * 0.001); // To avoid self intersection
        double dist_to_light = (light_pos - rec.p).length(); //As light is not at infinity

        hit_record shadow_rec;
        if (world->hit(shadow_ray_origin, l, 0.001, dist_to_light, shadow_rec)) {
            // We hit something! We are in shadow.
            return ambient; 
        }

        double diff = max(0.0, dot(rec.normal, l)); //To avoid dot prodect negetive
        vec3 diffuse = diff * rec.mat.color * light_color;

        double spec = pow(max(0.0, dot(rec.normal, h)), rec.mat.phong); // Exponent to make the reflection spot smaller
        vec3 specular = rec.mat.specular * spec * light_color;

        return ambient + diffuse + specular;
    }
    //If ray does not hit -> Surface behind something
    return vec3(0.2, 0.2, 0.2); // Dark grey background
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
    material red_plastic(vec3(0.7, 0.1, 0.1), 0.1, 0.9, 50.0);
    material green_matte(vec3(0.1, 0.8, 0.1), 0.5, 0.1, 1.0);
    material blue_metal(vec3(0.2, 0.2, 0.8), 0.2, 0.8, 20.0); 

    Surface* my_objects[3];

    // Object 0: The Target (Red Sphere)
    my_objects[0] = new Sphere(vec3(0, 0, -5), 1.0, red_plastic); 
    
    // Object 1: The Background (Green Triangle)
    my_objects[1] = new Triangle(vec3(-2,-2,-6), vec3(2,-2,-6), vec3(0,2,-6), green_matte);

    // Object 2: The Shadow Caster (Blue Sphere)
    // Placed between the light (10,10,10) and the Red Sphere
    my_objects[2] = new Sphere(vec3(1.5, 1.5, -2.5), 0.5, blue_metal);

    SurfaceList world(my_objects, 3);

    // 5. Render Loop (Outputting PPM format)
    // Header: P3 means ASCII color, then Width Height, then Max Color Value (255)
    cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    // Loop through every row (j) from Top to Bottom
    for (int j = image_height - 1; j >= 0; --j) {
        
        // Loop through every column (i) from Left to Right
        for (int i = 0; i < image_width; ++i) {
            
            // A. Calculate UV Coordinates (0.0 to 1.0)
            // This maps pixel (200, 200) to the center of the viewport (0.5, 0.5)
            double u = double(i) / (image_width - 1);
            double v = double(j) / (image_height - 1);

            // B. Create the Ray
            // Start at the eye (origin) and go through the specific spot on the viewport
            vec3 ray_direction = lower_left_corner + (u * horizontal) + (v * vertical) - origin;
            
            // C. Calculate Color (The Physics Engine)
            vec3 pixel_color = hit_color(origin, normalize(ray_direction), &world);

            // D. Write Color to Output
            // Convert 0.0-1.0 range to 0-255 integers
            cout << static_cast<int>(255.999 * pixel_color.x) << ' '
                 << static_cast<int>(255.999 * pixel_color.y) << ' '
                 << static_cast<int>(255.999 * pixel_color.z) << '\n';
        }
    }

    return 0;
}