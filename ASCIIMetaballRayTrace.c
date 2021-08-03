#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>     // For usleep()
#include <time.h>       // For clock() and CLOCKS_PER_SEC
#include <float.h>      // For DBL_MAX

#include <ncurses.h>    // compile with -lncurses
#include <math.h>       // compile with -lm

/* ASCIIMetaballRayTrace.c
 *
 * A basic metaball ray-tracer, rendering to the terminal
 *
 * Uses Ken Perlin's metaball function:
 *  Contribution from metaball i at point P = ( ((blobiness_i * (P - center_i)^2) / boundary_i^2) - blobiness_i )^exponent , where:
 *      - Blobiness affects the magnitude of a metaball's contribution; the center of metaball i contributes blobiness_i^exponent to the total.
 *      - The boundary of a metaball is the radius beyond which the metaball has no influence, i.e. the x-intercept of the contribution function.
 *      - The exponent must be a multiple of 2, and affects how steeply the curve falls. SPORE uses an exponent of 4, as does our current implementation.
 *
 *
 * TODO switch to proper sorted-intersection ray traversal; currently the bounding-sphere intersection computation is a means of pruning alone
 *
 * TODO screen Z-space still isn't standardized; the origin ought to be elsewhere, and the camera should do its own thing
 *
 * TODO shadows
 * December 2020
*/


// Vector-coordinate helpers
typedef struct vec3 {
    double x;
    double y;
    double z;
} vec3;


vec3 vec3Subtract(vec3 *a, vec3 *b){
    vec3 v = {a->x - b->x, a->y - b->y, a->z - b->z};
    return v;
}

double vec3DotProduct(vec3 *a, vec3 *b){
    return (a->x * b->x) + (a->y * b->y) + (a->z * b->z);
}

double vec3Length(vec3 *v){
    return sqrt(vec3DotProduct(v, v));
}

void vec3Normalize(vec3 *v){
    double length = vec3Length(v);
    v->x /= length;
    v->y /= length;
    v->z /= length;
}



// Luminance characters
#define LUMINANCE_CHARACTER_COUNT 12
const char *LUMINANCE_CHARACTERS = ".,:;-~=!*$#@";  // Nicked from the infamous donut.c


// Rendered world properties
#define SCREEN_WIDTH 100
#define SCREEN_HEIGHT 100
#define SCREEN_FROM_EYE_Z 200


vec3 LIGHT_SOURCE = {-100, 100, 0};  // Arbitrary
vec3 ORIGIN = {0, 0, 0};



// The metaballs
// TODO metaball algebra considerations: is the scale factor global or singular? how do the shaping func and threshold interact best?

#define METABALL_THRESHOLD 0.5  // Ranges from 0.0 (metaball will always fill bounding sphere) to (blobiness^exponent) or so, at which a metaball has no effect at its core

struct metaball {
    vec3 position;          // The position of the metaball in 3D space
    float blobiness;        // Scale factor for the metaball charge function; raises and lowers the y-intercept of the metaball charge quartic TODO rename
    float boundarySquared;  // Maximum radius of the metaball's area of effect, squared; moves the x-intercept of the metaball charge quartic (squared value)
};


// The metaball function (see file head)
double metaballEffectAtPoint(struct metaball *metaball, double radiusSquared){
    double effect = (metaball->blobiness * radiusSquared) / metaball->boundarySquared - metaball->blobiness;
    double tempThing = effect;  // TODO
    effect *= effect;   // ^2
   // effect *= effect;   // ^4
    effect *= tempThing;    // ^3   // TODO
    //return effect;
    return 0.0 - effect;    // TODO
}


#define METABALL_COUNT 5

struct metaball metaballs[METABALL_COUNT] = {
    {{0, 0, 300}, 1.0, 625},
    {{0, 0, 300}, 1.0, 10000},
    {{0, 0, 300}, 1.0, 2500},
    {{0, 0, 300}, 1.0, 1600},
    {{-50, -50, 350}, 1.0, 400}
};


// Function for computing the distance between two metaballs necessary to maintain a given radius at their intersection
double metaballDistanceFromRadius(struct metaball *a, struct metaball *b, double radius){
    // TODO doesn't work
    /*
    double aDist = sqrt((sqrt(METABALL_THRESHOLD/(2 * a->blobiness)) + 1) * a->boundarySquared - (radius * radius));
    double bDist = sqrt((sqrt(METABALL_THRESHOLD/(2 * b->blobiness)) + 1) * b->boundarySquared - (radius * radius));
    return aDist + bDist;
    */
    return sqrt(a->boundarySquared) *0.38 + sqrt(b->boundarySquared) *0.38; // func gives 0.5y at 0.4 x 
}

#define DESIRED_ARM_RADIUS 10   // Test value for the above




// Optimized geisswerks' algorithm entails sorting ray-metaball intercepts by distance; here, then, is an intercept
struct intercept {
    double length;
    int metaballIndex;
};

//...and a comparison function for said intercepts
int compareIntercepts(const void *a, const void *b){
    return ((struct intercept *)a)->length - ((struct intercept *)b)->length;
}


int main(void){
    // ncurses initialization sundries
    // More here: http://tldp.org/HOWTO/NCURSES-Programming-HOWTO/init.html
    initscr();           // Initialize the window
    noecho();            // Echo no keypresses
    curs_set(FALSE);     // Display no cursor
    keypad(stdscr,true); // Allow for arrow and function key
    nodelay(stdscr, true); // Set getch to non-blocking [STACK OVERFLOW sourced]

    int SCREEN_RESOLUTION_HEIGHT, SCREEN_RESOLUTION_WIDTH;
    getmaxyx(stdscr, SCREEN_RESOLUTION_HEIGHT, SCREEN_RESOLUTION_WIDTH);
    SCREEN_RESOLUTION_WIDTH /= 2;

    start_color();
    // Initialize colours
    init_pair(0, COLOR_RED, COLOR_BLACK);
    init_pair(1, COLOR_CYAN, COLOR_BLACK);
    init_pair(2, COLOR_GREEN, COLOR_BLACK);
    init_pair(3, COLOR_MAGENTA, COLOR_BLACK);
    init_pair(4, COLOR_YELLOW, COLOR_BLACK);




    // Compute the size of a pixel in world-space
    double PIXEL_WIDTH = ((double)SCREEN_WIDTH / (double)SCREEN_RESOLUTION_WIDTH) * ((double)SCREEN_RESOLUTION_WIDTH / (double)SCREEN_RESOLUTION_HEIGHT);
    double PIXEL_HEIGHT = (double)SCREEN_HEIGHT / (double)SCREEN_RESOLUTION_HEIGHT;

    double SCREEN_CENTER_X = ((double)SCREEN_WIDTH / 2.0) * ((double)SCREEN_RESOLUTION_WIDTH / (double)SCREEN_RESOLUTION_HEIGHT);
    double SCREEN_CENTER_Y = (double)SCREEN_HEIGHT / 2.0;


    int rotateBalls = 1;
    int rotateLight = 0;
    double prevTime = clock() / (double)CLOCKS_PER_SEC;
    double time = 0.0;
    double lightTime = 4.0;
    int done = 0;
    while(!done){

        clear();    // Clear screen

        // Iterate over every "pixel"
        for(int y = 0; y < SCREEN_RESOLUTION_HEIGHT; y++){
            for(int x = 0; x < SCREEN_RESOLUTION_WIDTH; x++){

                // STEP 1: Launch a ray from the eye into the screen, towards the given pixel
                // A ray is defined parametrically as an origin point (here, the eye; chosen to be (0, 0, 0))
                // and a normalized direction vector
                vec3 pixelPosition = {((double)x + 0.5) * PIXEL_WIDTH - SCREEN_CENTER_X, ((double)y + 0.5) * PIXEL_HEIGHT - SCREEN_CENTER_Y, SCREEN_FROM_EYE_Z};
                vec3 rayDirection = vec3Subtract(&pixelPosition, &ORIGIN);
                vec3Normalize(&rayDirection);
               


                // OPTIMIZED VARIANT
                // STEP 2: as per http://www.geisswerks.com/ryan/BLOBS/blobs.html, add metaballs whose bounding spheres intersect the ray to a list
                struct intercept intercepts[METABALL_COUNT * 2];
                int interceptCount = 0;
                
                // For free, we can establish ray distance boundaries (i.e. we don't check the ray from origin to screen, or beyond the farthest metaball)
                double maxT = 0.0;
                double minT = DBL_MAX;

                for(int i = 0; i < METABALL_COUNT; i++){
                    // Compute intercepts of ray and bounding sphere of metaball i
                    // We use a slightly modified form of the quadratic equation; see below
                    vec3 metaballCenterToOrigin = vec3Subtract(&ORIGIN, &(metaballs[i].position));
                    // a = 1.0; assuming ray direction is normalized
                    double b = 2.0 * vec3DotProduct(&metaballCenterToOrigin, &rayDirection);
                    double c = vec3DotProduct(&metaballCenterToOrigin, &metaballCenterToOrigin) - metaballs[i].boundarySquared;
                        
                    double discriminant = b * b - 4.0 * 1.0 * c;
                    if(discriminant > 0){
                        // Two roots; ray enters and exits bounding sphere; this is the only case we care about
                        // https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection says that
                        // this alternative version of the quadratic equation plays more nicely on machines, avoiding some float edge cases
                        double q = (b > 0) ? -0.5 * (b + sqrt(discriminant)) : -0.5 * (b - sqrt(discriminant));
                        
                        intercepts[interceptCount].length = q;
                        if (q > maxT) maxT = q;
                        intercepts[interceptCount].metaballIndex = i;
                        intercepts[interceptCount + 1].length = c / q;
                        if (c / q < minT) minT = c / q;
                        intercepts[interceptCount + 1].metaballIndex = i;
                        interceptCount += 2;
                    }

                }

                // We can use the above data as a primitive "pruner" to skip trace on rays that can't possibly intersect the blob
                // i.e. if there are intercepts, just run the unoptimized dirty version below
                if(interceptCount == 0){
                    //mvprintw(SCREEN_RESOLUTION_HEIGHT - y - 1, x*2, "XX");
                    continue;
                }

                /*
                // STEP 3: sort intercepts by distance
                // geisswerks source thinks mergesort is best (find out why)
                mergesort(intercepts, interceptCount, sizeof(struct intercept), &compareIntercepts);

                // STEP 4: iterate along ray, as in unoptimized version, but switching metaball contributions on and off based on intercept data
                // TODO
                */

                
                // UNOPTIMIZED DIRTY VERSION
                // Step along the ray in intervals
                for(double t = minT; t < maxT; t+= 0.5){    // Arbitrary step choices
                    // Compute point position of step
                    vec3 point;
                    point.x = ORIGIN.x + rayDirection.x * t;
                    point.y = ORIGIN.y + rayDirection.y * t;
                    point.z = ORIGIN.z + rayDirection.z * t;

                    // Sum metaball's charges
                    int closestMetaball;
                    double closestMetaballCharge = 0.0;
                    double charge = 0.0;
                    for(int i = 0; i < METABALL_COUNT; i++){   // for all metaballs
                        vec3 delta = vec3Subtract(&(metaballs[i].position), &point);
                        double rSquared = vec3DotProduct(&delta, &delta);
 
                        // Check bounding sphere
                        if(rSquared < metaballs[i].boundarySquared){
                            // if the point falls within the metaball's boundary, add the metaball's effect to the net effect at this point
                            double singleCharge = metaballEffectAtPoint(&metaballs[i], rSquared);
                            charge += singleCharge;
        
                            // Injected code for determining colour, etc.
                            if(singleCharge > closestMetaballCharge){
                                closestMetaballCharge = singleCharge;
                                closestMetaball = i;
                            }

                        }

                    }

                    // If charges are above threshold, we've entered the metaball blob
                    if(charge > METABALL_THRESHOLD){

                        // Compute normals based on all metaballs affecting this point
                        // Slightly redundant, as we only need to check metaballs in range and we should already know what those ones are
                        // TODO this is a step that would greatly benefit from the active-set optimization
                        vec3 normal = {0.0, 0.0, 0.0};
                        for(int i = 0; i < METABALL_COUNT; i++){   // for all metaballs
                            vec3 delta = vec3Subtract(&(metaballs[i].position), &point);
                            float rSquared = vec3DotProduct(&delta, &delta);
     
                            // Check bounding sphere (MAKE REDUNDANT!)
                            if(rSquared < metaballs[i].boundarySquared){
                                double singleCharge = metaballEffectAtPoint(&metaballs[i], rSquared);
                                // Compute normals (overdoing it here a bit; the charges should be done ONCE, at threshold above!)
                                vec3 singleNormal = vec3Subtract(&point, &(metaballs[i].position));
                                vec3Normalize(&singleNormal);
                                normal.x += singleNormal.x * singleCharge;
                                normal.y += singleNormal.y * singleCharge;
                                normal.z += singleNormal.z * singleCharge;
                            }

                        }
                        normal.x /= charge;
                        normal.y /= charge;
                        normal.z /= charge;



                        // Primitive normal lighting
                        vec3 vectorToLightSource = vec3Subtract(&LIGHT_SOURCE, &point);
                        vec3Normalize(&vectorToLightSource);
                        double luminance = vec3DotProduct(&normal, &vectorToLightSource);
                        double lumCharWid = (1.0 / LUMINANCE_CHARACTER_COUNT);
                        int interval;
                        if(luminance < 0){
                            interval = 0;
                        } else {
                            interval = (int)floor(luminance / lumCharWid);
                        }
                        // Color
                        attron(COLOR_PAIR(closestMetaball) | A_BOLD);
                        mvprintw(SCREEN_RESOLUTION_HEIGHT - y - 1, x*2, "%c%c", LUMINANCE_CHARACTERS[interval], LUMINANCE_CHARACTERS[interval]);
                        attroff(COLOR_PAIR(closestMetaball) | A_BOLD);
                        break;
                    }
                }
            }
        }

        refresh();  // Refresh the screen
        


        // Animate the metaballs
        double systemTime = clock() / (double)CLOCKS_PER_SEC;
        double deltaTime =  systemTime - prevTime;
        prevTime = systemTime;

        if(rotateBalls) time += deltaTime;
        if(rotateLight) lightTime += deltaTime;
        /*
        double zob;
        metaballs[0].position.x = -40.0; //0.0 - (100.0 - 100.0 * modf(time / 5.0, &zob));
        metaballs[0].position.y = 0.0;
        metaballs[0].position.z = 300.0;

        metaballs[1].position.x = 100.0 - 140.0 * modf(time / 20.0, &zob);
        metaballs[1].position.y = 0.0;
        metaballs[1].position.z = 300.0;
        */

        /*
        
        metaballs[0].position.z = 40.0 * sin(time * 1.3) + 300;
        metaballs[0].position.y = 40.0 * cos(time * 1.3) + 20;

        metaballs[1].position.z = 50.0 * sin(time) + 300;
        metaballs[1].position.x = 50.0 * cos(time) + 20.0 * sin(time * 0.5);
        metaballs[1].position.y = 10.0 * cos(time * 0.5);

        metaballs[2].position.z = 10.0 * sin(time * 1.8) + 300;
        metaballs[2].position.y = 10.0 * cos(time * 1.5);
        metaballs[2].position.x = 10.0 * sin(time * 1.2);

        double meta3radius = 50.0 + 10.0 * sin(time);
        metaballs[3].position.x = meta3radius * sin(time * 1.5);
        metaballs[3].position.y = meta3radius * cos(time * 1.5);
        metaballs[3].position.z = meta3radius * cos(time * 1.5) + 300;
        
        metaballs[4].blobiness = 50.0 * sin(time) + 60.0;
        */
        
        metaballs[0].position.x = 30.0 * sin(time * 1.83);
        metaballs[0].position.y = 30.0 * cos(time * 1.51);
        metaballs[0].position.z = 30.0 * cos(time * 1.22) + 400.0;

        metaballs[1].position.x = metaballs[0].position.x + metaballDistanceFromRadius(metaballs+0, metaballs+1, DESIRED_ARM_RADIUS) * sin(time * 0.8);
        metaballs[1].position.y = metaballs[0].position.y + metaballDistanceFromRadius(metaballs+0, metaballs+1, DESIRED_ARM_RADIUS) * cos(time * 0.8);
        metaballs[1].position.z = metaballs[0].position.z +  metaballDistanceFromRadius(metaballs+0, metaballs+1, DESIRED_ARM_RADIUS)* sin(time * 0.8);

        metaballs[2].position.x = metaballs[1].position.x + metaballDistanceFromRadius(metaballs+1, metaballs+2, DESIRED_ARM_RADIUS) * sin(time * 2.2);
        metaballs[2].position.y = metaballs[1].position.y + metaballDistanceFromRadius(metaballs+1, metaballs+2, DESIRED_ARM_RADIUS) * cos(time * 2.2);
        metaballs[2].position.z = metaballs[1].position.z + metaballDistanceFromRadius(metaballs+1, metaballs+2, DESIRED_ARM_RADIUS) * sin(time * 2.2);

        metaballs[3].position.x = metaballs[2].position.x + metaballDistanceFromRadius(metaballs+2, metaballs+3, DESIRED_ARM_RADIUS) * sin(time * 0.2);
        metaballs[3].position.y = metaballs[2].position.y + metaballDistanceFromRadius(metaballs+2, metaballs+3, DESIRED_ARM_RADIUS) * cos(time * 0.2);
        metaballs[3].position.z = metaballs[2].position.z + metaballDistanceFromRadius(metaballs+2, metaballs+3, DESIRED_ARM_RADIUS) * sin(time * 0.2);

        metaballs[4].position.x = metaballs[3].position.x + metaballDistanceFromRadius(metaballs+3, metaballs+4, DESIRED_ARM_RADIUS) * sin(time * 1.8);
        metaballs[4].position.y = metaballs[3].position.y + metaballDistanceFromRadius(metaballs+3, metaballs+4, DESIRED_ARM_RADIUS) * cos(time * 1.8);
        metaballs[4].position.z = metaballs[3].position.z + metaballDistanceFromRadius(metaballs+3, metaballs+4, DESIRED_ARM_RADIUS) * sin(time * 1.8);


        // Animate the light
        LIGHT_SOURCE.x = 100.0 * cos(lightTime);
        LIGHT_SOURCE.y = 100.0;
        LIGHT_SOURCE.z = 100.0 * sin(lightTime) + 300;
        //LIGHT_SOURCE = metaballs[3].position;
        
        // Check for user input
        switch(getch()){
            case ' ':
                rotateBalls = !rotateBalls;
            break;
            case 'l':
                rotateLight = !rotateLight;
            break;
            case 'x':
                done = 1;
            break;
        }

        // Don't overdo it
        usleep(1000);

    }
    endwin();
}
