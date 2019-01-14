#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

// Macros --------------------------------------------------------------------------------------------------------------
#define REAL_MSL_TRAJECTORY

// Data Types ----------------------------------------------------------------------------------------------------------
// Type shortcut
typedef unsigned long long int  super_int;

// High-precision storage container
typedef sf::Vector3<long double> Vector3d;

// Celestial body list
enum Bodies {
    SUN,
    //MERCURY,
    VENUS,
    EARTH,
    MOON,
    MARS,
    JUPITER,
    //SATURN,
    //URANUS,
    //NEPTUNE,
    NUM_BODIES
};

// Struct for initial celestial body data
struct Init_Body_Data {
    double mass;            // in kg
    Vector3d position;      // in m
    Vector3d velocity;      // in m/s
    float radius;           // pixel radius
    sf::Color colour;       // display colour
    std::string name;       // display name
};

// Constants -----------------------------------------------------------------------------------------------------------
// General Constants
const long double   MINUTE              = 60;                       // Seconds per minute
const long double   HOUR                = 3600;                     // Seconds per hour
const long double   DAY                 = 86400;                    // Seconds per day
const long double   YEAR                = 31536000;                 // Seconds per year

// Calculation Constants
const long double   GRAV_CONSTANT       = 6.67408e-11;              // Gravitational constant in m^3 kg^-1 s^-2
const long double   POSITION_SCALE      = 1e-7;                     // 1 pixel per 10,000,000 m
const long double   TIME_STEP           = 1 * HOUR;                 // Time interval for calculations
const long double   SIM_DURATION        = 3 * YEAR;                 // Sim duration in sim time, 1 year simulation
const super_int     NUM_TIME_STEPS      = SIM_DURATION / TIME_STEP; // Number of time steps in the simulation

// Launch Velocity-Specific Constants
const long double   LAUNCH_SPEED        = 36200;                    // Speed obtainable by spaceship in m/s
const int           NUM_BRACKETING_PTS  = 10;                       // Number of bracketing points in each x/y direction
const long double   ARRIVAL_DISTANCE    = 1e8;                      // Maximum acceptable arrival distance in m

// Simulation Constants
const int           WINDOW_X            = 1200;                     // Window size
const int           WINDOW_Y            = 600;
const float         HUD_X               = 300;                      // HUD size
const float         HUD_Y               = 150;
const float         DEFAULT_ZOOM        = 50;                       // Default zoom
const long double   DEFAULT_TIME_SCALE  = 60 * TIME_STEP;           // Default time scale, 60 hours per second
const float         DEFAULT_PAN_SPEED   = 50;                       // Default pan speed
const float         ZOOM_SCROLL_FACTOR  = 0.1;                      // Zoom speed with scrolling
const float         SIZE_SCROLL_FACTOR  = 0.75;                     // Resize speed from zooming
const long double   TIME_SCROLL_FACTOR  = 0.1;                      // Time scale change speed with scrolling

// Celestial Data Constants
Init_Body_Data BODY_DATA[NUM_BODIES] = {
    [SUN] = {
        .mass       = 1.988544e30,
        .position   = Vector3d(-5.067082710407100e8, -1.738890230441533e8, 2.998128794486838e5),
        .velocity   = Vector3d(6.829281815579370e0, -8.508243568302674e0, -1.322951483997670e-1),
        .radius     = 300,
        .colour     = sf::Color::Yellow,
        .name       = "SUN",
    },

//    [MERCURY] = {
//        .mass       = 3.302e23,
//        .position   = Vector3d(4.385121966632131e10, 2.106068877817014e10, -2.334757852275234e9),
//        .velocity   = Vector3d(-3.051639850372996e4, 4.608448778852690e4, 6.566560938576007e3),
//        .radius     = 30,
//        .colour     = sf::Color(230, 210, 205),
//        .name       = "MERCURY",
//    },

    [VENUS] = {
        .mass       = 48.685e23,
        .position   = Vector3d(6.455211776637389e10, -8.740761317312035e10, -4.949628397207707e9),
        .velocity   = Vector3d(2.784762776612149e4, 2.080557195687820e4, -1.321747469962308e3),
        .radius     = 85,
        .colour     = sf::Color(255, 160, 75),
        .name       = "VENUS",
    },

    [EARTH] = {
        .mass       = 5.97219e24,
        .position   = Vector3d(6.451911861931717e10, 1.323830938644331e11, -3.743509191088378e6),
        .velocity   = Vector3d(-2.723506168441186e4, 1.300170564145261e4, -1.650464599978640e0),
        .radius     = 100,
        .colour     = sf::Color::Blue,
        .name       = "EARTH",
    },

    [MOON] = {
        .mass       = 734.9e20,
        .position   = Vector3d(6.447516810277499e10, 1.320198456044836e11, 1.573990353137255e6),
        .velocity   = Vector3d(-2.617767382646868e4, 1.282353975402548e4, 9.695532302813259e1),
        .radius     = 15,
        .colour     = sf::Color(200, 200, 200),
        .name       = "MOON",
    },

    [MARS] = {
        .mass       = 6.4185e23,
        .position   = Vector3d(-1.232095936109405e11, 2.111037367449816e11, 7.440075911190972e9),
        .velocity   = Vector3d(-2.002706327037775e4, -1.011594824595273e4, 2.799958469180357e2),
        .radius     = 60,
        .colour     = sf::Color::Red,
        .name       = "MARS",
    },

    [JUPITER] = {
        .mass       = 1898.13e24,
        .position   = Vector3d(5.871490243900723e11, 4.543647477078475e11, -1.503783743410757e10),
        .velocity   = Vector3d(-8.151911318827203e3, 1.095839062614485e4, 1.369343353961496e2),
        .radius     = 150,
        .colour     = sf::Color(255, 185, 115),
        .name       = "JUPITER",
    },

//    [SATURN] = {
//        .mass       = 5.68319e26,
//        .position   = Vector3d(-1.347394185707304e12, -5.292083712121876e11, 6.282302197861078e10),
//        .velocity   = Vector3d(3.010364375775313e3, -9.013422046378633e3, 3.756415821239356e1),
//        .radius     = 130,
//        .colour     = sf::Color(255, 200, 50),
//        .name       = "SATURN",
//    },
//
//    [URANUS] = {
//        .mass       = 86.8103e24,
//        .position   = Vector3d(2.998541955162921e12, 1.605015560992119e11, -3.825181389328766e10),
//        .velocity   = Vector3d(-4.138273607890305e2, 6.482743430644901e3, 2.935685778322350e1),
//        .radius     = 130,
//        .colour     = sf::Color(100, 235, 160),
//        .name       = "URANUS",
//    },
//
//    [NEPTUNE] = {
//        .mass       = 102.41e24,
//        .position   = Vector3d(3.883498525853999e12, -2.249388560079292e12, -3.825181389328766e10),
//        .velocity   = Vector3d(2.688803868386490e3, 4.735311073379322e3, 2.935685778322350e1),
//        .radius     = 130,
//        .colour     = sf::Color::Blue,
//        .name       = "NEPTUNE",
//    },
};

Init_Body_Data SPACE_PROBE_DATA = {
    .mass       = 3893,
    .position   = Vector3d(6.452460147257443e10,  1.323861458716573e11, -1.093350118016452e7),
    #ifdef REAL_MSL_TRAJECTORY
    .velocity   = Vector3d(-2.928870616986328E+04, 2.046002445594545E+04, -5.867797836085259E+03),
    #else
    .velocity   = Vector3d(0, 0, 0),
    #endif
    .radius     = 100,
    .colour     = sf::Color::Magenta,
    .name       = "SPACE PROBE",
};

// Globals -------------------------------------------------------------------------------------------------------------
sf::Font font;

// Functions -----------------------------------------------------------------------------------------------------------
// Calculate distance between two (x,y,z) coordinates (vectors)
long double vector_distance(Vector3d & vec_1, Vector3d & vec_2) {
    return (long double)sqrt((vec_1.x - vec_2.x)*(vec_1.x - vec_2.x) +
                             (vec_1.y - vec_2.y)*(vec_1.y - vec_2.y) +
                             (vec_1.z - vec_2.z)*(vec_1.z - vec_2.z));
}

// Classes -------------------------------------------------------------------------------------------------------------
// Celestial body class
class Body : public sf::Drawable, public sf::Transformable {
public:
    // Data constructor
    Body(Init_Body_Data *init_data) : real_position(init_data->position), velocity(init_data->velocity),
                                      shape(init_data->radius) {
        mass = init_data->mass;
        shape.setOrigin(init_data->radius, init_data->radius);
        shape.setFillColor(init_data->colour);
        name.setFont(font);
        name.setString(init_data->name);
        name.setOrigin(init_data->radius, init_data->radius);
        name.setFillColor(sf::Color::White);
        name.setCharacterSize(500);
        name.setPosition(init_data->radius * 1.5, init_data->radius * -1.5);
    }

    // Calculate acceleration caused by another body
    void Calculate_Acc(Body * other) {
        acceleration += Vector3d(other->real_position - real_position) * GRAV_CONSTANT * other->mass /
                        (long double)pow(vector_distance(real_position, other->real_position), 3);
    }

    // Calculate velocity and position
    void Calculate_Vel_Pos(long double time_elapsed) {
        velocity += acceleration * time_elapsed;
        real_position += velocity * time_elapsed;
        acceleration = Vector3d(0, 0, 0);
    }

    // Manually set the real position
    void SetRealPosition(Vector3d & new_position) {
        real_position = new_position;
    }

    // Manually set the velocity
    void SetVelocity(Vector3d & new_velocity) {
        velocity = new_velocity;
    }

    // Fetch the real position
    Vector3d GetRealPosition() const {
        return real_position;
    }

    // Fetch the velocity
    Vector3d GetVelocity() const {
        return velocity;
    }

private:
    // Render the shape
    void draw(sf::RenderTarget & target, sf::RenderStates states) const {
        states.transform *= this->getTransform();
        target.draw(shape, states);
        target.draw(name, states);
    }

    // Private members
    Vector3d real_position, velocity, acceleration;
    sf::CircleShape shape;
    sf::Text name;
    long double mass;
};

// Main Function -------------------------------------------------------------------------------------------------------
int main() {
    // Parameter display -----------------------------------
    std::cout << "Welcome to the Mars mission trajectory calculator.\n"
              << "-------------------------\n\n"
              << "Launch date and time: 2011-11-26, 15:52 UTC.\n"
              << "-------------------------\n\n"
              << "Simulation paramters:"
              << "\n\tSimulation duration: " << SIM_DURATION / DAY << " days"
              << "\n\tTime step: " << TIME_STEP / DAY << " days"
              << "\n\tNumber of bracketing points every iteration: " << NUM_BRACKETING_PTS << " divisions"
              << "\n\tLaunch speed: " << LAUNCH_SPEED << " m/s"
              << "\n\tAcceptable arrival distance: " << ARRIVAL_DISTANCE / 1000 << " km"
              << "\n-------------------------\n\n";

    // Body Initialization----------------------------------
    std::cout << "Initializing bodies...\n";

    Body *body[NUM_BODIES];
    for (int i = 0; i < NUM_BODIES; i++) {
        body[i] = new Body(&BODY_DATA[i]);
    }
    Body Space_Probe(&SPACE_PROBE_DATA);

    // Celestial Body Calculations -------------------------
    std::cout << "Calculating positions of celestial bodies...\n";

    // Calculate each body's position at every time interval and store the data
    static Vector3d position[NUM_TIME_STEPS][NUM_BODIES]; // Very large array
    for (super_int time_index = 0; time_index < NUM_TIME_STEPS; time_index++) {
        for (int i = 0; i < NUM_BODIES; i++) {
            for (int j = 0; j < NUM_BODIES; j++) {
                if (i != j) {
                    body[i]->Calculate_Acc(body[j]);
                }
            }
        }

        for (int i = 0; i < NUM_BODIES; i++) {
            body[i]->Calculate_Vel_Pos(TIME_STEP);
            position[time_index][i] = body[i]->GetRealPosition();
        }
    }

    // Spaceship Trajectory Calculations -------------------
    // Pointer to the calculated trajectory
    Vector3d *space_probe_trajectory;

    #ifdef REAL_MSL_TRAJECTORY //---------------------------------------------------------------------------------------
    // If simulated with the real initial velocity:

    std::cout << "Calculating the real MSL trajectory...\n";

    // Storage array for the position
    static Vector3d space_probe_position[NUM_TIME_STEPS];
    space_probe_position[0] = SPACE_PROBE_DATA.position;

    // Set the velocity
    space_probe_trajectory = &space_probe_position[0];

    // Iterate through each time index to calculate the trajectory
    for (super_int time_index = 0; time_index < NUM_TIME_STEPS - 1; time_index++) {
        // Get each body's position at the time index, and the probe's acceleration due to that body
        for (int c = 0; c < NUM_BODIES; c++) {
            body[c]->SetRealPosition(position[time_index][c]);
            Space_Probe.Calculate_Acc(body[c]);
        }

        // Calculate and store the probe position data
        Space_Probe.Calculate_Vel_Pos(TIME_STEP);
        space_probe_position[time_index + 1] = Space_Probe.GetRealPosition();
    }

    std::cout << "Complete.\n"
              << "-------------------------\n\n";

    #else // -----------------------------------------------------------------------------------------------------------
    // If calculating the initial velocity:

    std::cout << "Calculating and comparing spaceship trajectories...\n";
    // Storage arrays for position and velocity

    static Vector3d space_probe_position[NUM_BRACKETING_PTS][NUM_BRACKETING_PTS][2][NUM_TIME_STEPS];
    static Vector3d space_probe_velocity[NUM_BRACKETING_PTS][NUM_BRACKETING_PTS][2][NUM_TIME_STEPS];

    // Trajectory variables
    sf::Vector3i best_launch_params(0, 0, 0);
    float time_of_interception;
    long double min_mars_distance = 1e15, temp_mars_distance;

    // Bracketing variables
    Vector3d bracketing_center(0, 0, 0);
    long double bracketing_bound = LAUNCH_SPEED;
    int bracketing_iteration = 0;
    long double temp_vel_x, temp_vel_y;

    // Search for a launch vector that gives an acceptable distance to Mars
    while (min_mars_distance > ARRIVAL_DISTANCE) {
        // Generate an initial x,y,z launch velocity vector with the launch speed at each bracketing point
        for (int i = 0; i < NUM_BRACKETING_PTS; i++) {
            // Iterate through each x bracketing point within bracketing_bound from bracketing_center
            temp_vel_x = bracketing_center.x + bracketing_bound * (-1 + 2 * (float)i / (NUM_BRACKETING_PTS - 1));
            for (int j = 0; j < NUM_BRACKETING_PTS; j++) {
                // Iterate through each y bracketing point within bracketing_bound from bracketing_center
                temp_vel_y = bracketing_center.y + bracketing_bound * (-1 + 2 * (float)j / (NUM_BRACKETING_PTS - 1));
                for (int k = 0; k < 2; k++) {
                    space_probe_position[i][j][k][0] = SPACE_PROBE_DATA.position;
                    space_probe_velocity[i][j][k][0].x = temp_vel_x;
                    space_probe_velocity[i][j][k][0].y = temp_vel_y;
                    // Calculate the z-coordinate for this vector
                    space_probe_velocity[i][j][k][0].z = -k * sqrt(LAUNCH_SPEED * LAUNCH_SPEED -
                                                              space_probe_velocity[i][j][k][0].x *
                                                              space_probe_velocity[i][j][k][0].x -
                                                              space_probe_velocity[i][j][k][0].y *
                                                              space_probe_velocity[i][j][k][0].y);
                }
            }
        }

        // Iterate through each time index to calculate the trajectory
        for (super_int time_index = 0; time_index < NUM_TIME_STEPS - 1; time_index++) {
            // Get each body's position at the time index
            for (int c = 0; c < NUM_BODIES; c++) {
                body[c]->SetRealPosition(position[time_index][c]);
            }

            // Apply the dynamic forces to the spacecraft's trajectory from each launch vector
            for (int i = 0; i < NUM_BRACKETING_PTS; i++) {
                for (int j = 0; j < NUM_BRACKETING_PTS; j++) {
                    for (int k = 0; k < 2; k++) {
                        // Discard the trajectory if it had an invalid launch vector
                        if (sqrt(space_probe_velocity[i][j][k][0].x * space_probe_velocity[i][j][k][0].x +
                                 space_probe_velocity[i][j][k][0].y * space_probe_velocity[i][j][k][0].y) <=
                            LAUNCH_SPEED) {
                            // Load the position and velocity for the trajectory at the time index
                            Space_Probe.SetRealPosition(space_probe_position[i][j][k][time_index]);
                            Space_Probe.SetVelocity(space_probe_velocity[i][j][k][time_index]);

                            // Apply acceleration, velocity, and position changes
                            for (int c = 0; c < NUM_BODIES; c++) {
                                Space_Probe.Calculate_Acc(body[c]);
                            }
                            Space_Probe.Calculate_Vel_Pos(TIME_STEP);

                            // Save the position and velocity data
                            space_probe_position[i][j][k][time_index + 1] = Space_Probe.GetRealPosition();
                            space_probe_velocity[i][j][k][time_index + 1] = Space_Probe.GetVelocity();

                            // Calculate and compare the distance to Mars
                            temp_mars_distance = vector_distance(space_probe_position[i][j][k][time_index + 1],
                                                                 position[time_index + 1][MARS]);
                            if (temp_mars_distance < min_mars_distance) {
                                min_mars_distance = temp_mars_distance;
                                best_launch_params = sf::Vector3i(i, j, k);
                                time_of_interception = time_index * TIME_STEP;
                            }
                        }
                    }
                }
            }
        }

        // Output the bracketing data
        bracketing_iteration++;
        std::cout << "\tIteration " << bracketing_iteration << ": " << min_mars_distance / 1000 << " km\n";

        // Update the bracketing center and bound
        bracketing_center = space_probe_velocity[best_launch_params.x][best_launch_params.y][best_launch_params.z][0];
        bracketing_bound /= NUM_BRACKETING_PTS - 1;
    }

    // Set the trajectory pointer to the calculated position
    space_probe_trajectory = &space_probe_position[best_launch_params.x][best_launch_params.y][best_launch_params.z][0];

    // Output the simulation data
    std::cout << "Complete.\n"
              << "-------------------------\n\n"
              << "Best launch vector, in m/s:"
              << "\n\tx: " << bracketing_center.x
              << "\n\ty: " << bracketing_center.y
              << "\n\tz: " << bracketing_center.z
              << "\nIt will pass within " << min_mars_distance / 1000 << " km of Mars, after "
              << time_of_interception / DAY << " days of travel.\n"
              << "-------------------------\n\n";

    #endif // ----------------------------------------------------------------------------------------------------------

    // Simulation Initialization ---------------------------
    std::cout << "Preparing simulation...\n";

    // Rendering window
    sf::RenderWindow window(sf::VideoMode(WINDOW_X, WINDOW_Y), "N-Body Simulation", sf::Style::Close);

    // Font set-up
    font.loadFromFile("OCRAEXT.ttf");

    // Main view
    sf::View main_view(sf::Vector2f(0, 0), sf::Vector2f(WINDOW_X, WINDOW_Y));
    main_view.zoom(DEFAULT_ZOOM);
    float pan_speed = DEFAULT_PAN_SPEED;

    // HUD view
    sf::View hud_view(sf::FloatRect(0, 0, HUD_X, HUD_Y));
    hud_view.setViewport(sf::FloatRect((WINDOW_X - HUD_X) / WINDOW_X, (WINDOW_Y - HUD_Y) / WINDOW_Y,
                                       HUD_X / WINDOW_X, HUD_Y / WINDOW_Y));

    // HUD time display
    sf::Text time_display;
    time_display.setFont(font);
    time_display.setPosition(0, 0);
    time_display.setFillColor(sf::Color::White);
    time_display.setCharacterSize(30);
    std::stringstream time_display_SS;
    time_display_SS << std::setprecision(4);

    // Clock
    sf::Clock clock;
    long double time_elapsed = 0;
    super_int time_elapsed_index = 0;
    long double time_scale = DEFAULT_TIME_SCALE;

    // Mouse
    sf::Vector2i mouse_prev_pos;
    bool mouse_prev_pressed;

    // Other flags and counters
    bool sim_running = false;

    std::cout << "Simulating.\n"
              << "-------------------------\n\n";

    // Main Loop -------------------------------------------
    while (window.isOpen()) {
        // Event Handling ----------------------------------
        sf::Event event;
        while (window.pollEvent(event)) {
            // Window closed
            if (event.type == sf::Event::Closed) {
                window.close();
            }

            // Pause/unpause simulation
            if (event.type == sf::Event::MouseButtonPressed) {
                if (event.mouseButton.button == sf::Mouse::Right) {
                    // Pause/unpause simulation
                    sf::Time pause_cycle = clock.restart();
                    if (sim_running) {
                        time_elapsed += pause_cycle.asSeconds() * time_scale;
                    }
                    sim_running = !sim_running;
                }
                else if (event.mouseButton.button == sf::Mouse::Middle) {
                    // Reverse time
                    time_scale *= -1;
                }
            }

            if (event.type == sf::Event::MouseWheelScrolled) {
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::LControl)) {
                    // Time control
                    time_scale *= 1 + event.mouseWheelScroll.delta * TIME_SCROLL_FACTOR;
                }
                else {
                    // Zoom in/out
                    float curr_zoom_factor = event.mouseWheelScroll.delta * ZOOM_SCROLL_FACTOR;
                    // Change the zoom and pan speeds
                    main_view.zoom(1 + curr_zoom_factor);
                    pan_speed *= 1 + curr_zoom_factor;

                    // Change the body display sizes for visibility
                    for (int i = 0; i < NUM_BODIES; i++) {
                        body[i]->scale(1 + curr_zoom_factor * SIZE_SCROLL_FACTOR,
                                       1 + curr_zoom_factor * SIZE_SCROLL_FACTOR);
                    }

                    Space_Probe.scale(1 + curr_zoom_factor * SIZE_SCROLL_FACTOR,
                                 1 + curr_zoom_factor * SIZE_SCROLL_FACTOR);
                }
            }
        }

        // Detailed Mouse Handling -------------------------
        // Panning (drag)
        if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
            if (mouse_prev_pressed) {
                // Pan through the main view
                main_view.move(static_cast<sf::Vector2f>(sf::Mouse::getPosition(window) - mouse_prev_pos) * -pan_speed);
            }
            else {
                // Start panning
                window.setMouseCursorVisible(false);
                mouse_prev_pressed = true;
            }

            mouse_prev_pos = sf::Mouse::getPosition(window);
        }
        else if (mouse_prev_pressed) {
            // Stop panning
            window.setMouseCursorVisible(true);
            mouse_prev_pressed = false;
        }

        // Simulation Logic --------------------------------
        if (sim_running) {
            // Calculate the simulation time
            time_elapsed += clock.restart().asSeconds() * time_scale;
            time_elapsed_index = time_elapsed / TIME_STEP;

            // Ensure the simulation is within the constraints
            if (time_elapsed < 0) {
                time_elapsed = 0;
                time_scale = fabs(time_scale);
                sim_running = false;
            }
            else if (time_elapsed > SIM_DURATION) {
                time_elapsed = SIM_DURATION - TIME_STEP;
                sim_running = false;
            }

            #ifndef REAL_MSL_TRAJECTORY
            // Pause the simulation if intercepted Mars
            if ((super_int)(time_of_interception / TIME_STEP) == time_elapsed_index) {
                sim_running = false;
            }
            #endif

            // Fetch each body's position at the simulation time index
            for (int i = 0; i < NUM_BODIES; i++) {
                body[i]->setPosition(position[time_elapsed_index][i].x * POSITION_SCALE,
                                     -position[time_elapsed_index][i].y * POSITION_SCALE);
            }

            Space_Probe.setPosition(space_probe_trajectory[time_elapsed_index].x * POSITION_SCALE,
                                    -space_probe_trajectory[time_elapsed_index].y * POSITION_SCALE);
        }

        // HUD Update --------------------------------------
        // Feed the simulation time into a string
        time_display_SS.str(std::string());
        time_display_SS << time_elapsed / YEAR << "\nyears\nelapsed";
        #ifndef REAL_MSL_TRAJECTORY
        if ((super_int)(time_of_interception / TIME_STEP) == time_elapsed_index) {
            time_display_SS << "\nInterception";
        }
        #endif
        time_display.setString(time_display_SS.str());

        // Graphics Rendering ------------------------------
        window.clear();

        // Render into the main view
        window.setView(main_view);
        for (int i = 0; i < NUM_BODIES; i++) {
            window.draw(*body[i]);
        }
        window.draw(Space_Probe);

        // Render into the HUD
        window.setView(hud_view);
        window.draw(time_display);

        window.display();
    }

    // De-initialization -----------------------------------
    std::cout << "Simulation complete.\n"
              << "De-initializing...\n";
    for (int i = 0; i < NUM_BODIES; i++) {
        delete body[i];
    }

    std::cout << "Finished!\n"
              << "-------------------------\n\n";
    std::cin.get();
    return 0;
}

