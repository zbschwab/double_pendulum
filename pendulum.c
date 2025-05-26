#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <SDL.h>
#include <SDL_image.h>
#include <SDL_render.h>

#include "phys_math.h"

#define MAXTIME 60  // in seconds
#define T_STEP 0.01  // timestep dt in seconds
#define PI 3.1415927
#define G 9.8

#define CHARLIMIT 10
#define SCREEN_WIDTH 640
#define SCREEN_HEIGHT 480
#define FRAME_RATE 60
#define SCALE_CONST 100  // make calculated coords a reasonable size in pixels
#define PIVOT_SIZE 10

// credit: SDL setup from Charlie's lab @ https://curtsinger.cs.grinnell.edu/teaching/2025S/CSC161/readings/sdl2/              
//              and @ https://git.cs.grinnell.edu/csc161-s25/SDL2-demo/src/commit/a6ce6b007d9409db9c3b4e882cc95f6690ea3b89/demo3.c
//         SDL wiki referenced for writing animation in general @ https://wiki.libsdl.org/SDL2/CategoryAPIFunction
//         example code referenced for SDL rendering images @ https://gist.github.com/armornick/3434362
//         double pendulum derivation (my own work) @ https://www.overleaf.com/read/pnfwfhzjxhsd#2c25da
//              (checked with lagrangian double pendulum results @ dassencio.org/33)
//         Runge-Kutta 4th-order 2nd-order-ODE solver written with help from wikipedia and a fortran screenshot
//              https://math.stackexchange.com/q/3095319
//              https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

// mechanics precision conventions: rad needs 6 digits, kg and meter need 3

// note: this program leaks 184 bytes from 1 object. i'm pretty sure this isn't me, i
//       think the low-level library 'lidbus' SDL uses for some linux architectures leaks.

double* get_init_conds() {
    char user_input[CHARLIMIT];
    int input_len;
    double* init_conds = (double*)malloc(6*sizeof(double));
    char* ic_prompts[] = {"Mass of 1st pendulum (0<m<10 kg): ",
                          "Mass of 2nd pendulum (0<m<10 kg): ",
                          "Length of 1st pendulum (10<l<100 cm): ",
                          "Length of 2nd pendulum (10<l<100 cm): ",
                          "Initial angle of 1st pendulum (0<a<180 deg)",
                          "Initial angle of 2nd pendulum (0<a<180 deg)"};

    printf("Welcome to the double pendulum simulator. Specify your initial conditions:\n");

    for (int i = 0; i < 6; i++) {
        printf("%s\n", ic_prompts[i]);
        if (fgets(user_input, CHARLIMIT, stdin) == NULL) {
            fprintf(stderr, "Failed to read input from user.\n");
            exit(EXIT_FAILURE);
        }
        input_len = strlen(user_input);
        // check if user input is valid length
        if (input_len > 6) {
            printf("Invalid input: too long. Please try again.\n");
            return get_init_conds();
        } else if (input_len == 0) {
            printf("Invalid input: enter a number. Please try again.\n");
            return get_init_conds();
        }
        char* endptr;
        double ic = strtod(user_input, &endptr);  // convert input to double
        if ((ic == 0) && (user_input[0] != '0')) {
            printf("Invalid input: non-number character. Please try again.\n");
            return get_init_conds();
        }
        init_conds[i] = ic;
    }

    // check number range
    if ((init_conds[M1] <= 0) | (init_conds[M1] > 10) | (init_conds[M2] <= 0) | (init_conds[M2] > 10)) {
        printf("Invalid input: mass out of range, Please try again.\n");
        return get_init_conds();
    } else if ((init_conds[L1] < 10) | (init_conds[L1] > 100) | (init_conds[L2] < 10) | (init_conds[L2] > 100)) {
        printf("Invalid input: length out of range, Please try again.\n");
        return get_init_conds();
    } else if ((init_conds[T1] < 0) | (init_conds[T1] > 180) | (init_conds[T2] < 0) | (init_conds[T2] > 180)) {
        printf("Invalid input: angle out of range, Please try again.\n");
        return get_init_conds();
    }

    return init_conds;
}

int main(void) {

    double* init_conds = get_init_conds();

    // put user-specified constant values in struct
    constants_t *constants = (constants_t*)malloc(sizeof(constants_t));
    constants->m_1 = init_conds[M1];
    constants->m_2 = init_conds[M2];
    constants->l_1 = init_conds[L1] / 100;  // convert cm to m
    constants->l_2 = init_conds[L2] / 100;

    // initialize system state with user-specified initial angles
    state_t *state = (state_t*)malloc(sizeof(state_t));
    state->theta_1 = deg_to_rad(init_conds[T1]);
    state->theta_2 = deg_to_rad(init_conds[T2]);
    state->omega_1 = 0;
    state->omega_2 = 0;

    // total number of timesteps for the defined duration
    int total_steps = (int)round(MAXTIME/T_STEP);

    // declare array to hold angle values over time in 
    double theta_1_arr[total_steps];
    double theta_2_arr[total_steps];
    memset(theta_1_arr, 0, total_steps*sizeof(double));
    memset(theta_2_arr, 0, total_steps*sizeof(double));

    // fill 1st slot in array with initial angles
    theta_1_arr[0] = init_conds[T1];
    theta_2_arr[0] = init_conds[T2];

    double t = 0;

    // call RK4
    for (int i = 0; i < total_steps; i++) {
        runge_kutta_4(state, constants, t);
        theta_1_arr[i] = state->theta_1;
        theta_2_arr[i] = state->theta_2;
        t += T_STEP;
    }

    // arrays to hold pendulum angles in cartesian coordinates 
    double x_1_arr[total_steps];
    double x_2_arr[total_steps];
    double y_1_arr[total_steps];
    double y_2_arr[total_steps];
    memset(x_1_arr, 0, total_steps*sizeof(double));
    memset(x_2_arr, 0, total_steps*sizeof(double));
    memset(y_1_arr, 0, total_steps*sizeof(double));
    memset(y_2_arr, 0, total_steps*sizeof(double));

    cart_t *cart = (cart_t*)malloc(sizeof(cart_t));

    // convert pendulum angles in polar to cartesian
    for (int i = 0; i < total_steps; i++) {
        polar_to_car(theta_1_arr[i], SCALE_CONST*constants->l_1, cart);
        x_1_arr[i] = cart->x;
        y_1_arr[i] = cart->y;
        polar_to_car(theta_2_arr[i], SCALE_CONST*constants->l_2, cart);
        x_2_arr[i] = cart->x + x_1_arr[i];
        y_2_arr[i] = cart->y + y_1_arr[i];
    }

    // simulate double pendulum motion with SDL
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        // failed to initialize SDL, log error, exit program
        fprintf(stderr, "Failed to initialize SDL: %s\n", SDL_GetError());
        exit(EXIT_FAILURE);
    }

    // create SDL window (centered on screen, 640x480)
    SDL_Window* window = SDL_CreateWindow("Double Pendulum Simulator",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        SCREEN_WIDTH, SCREEN_HEIGHT, 
        SDL_WINDOW_SHOWN);
    if (window == NULL) {
        // failed to initialize window, log error, shut down SDL, exit
        fprintf(stderr, "Failed to create SDL window: %s\n", SDL_GetError());
        SDL_Quit();
        exit(EXIT_FAILURE);
    }

    // create a renderer for the window
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (renderer == NULL) {
        // detect error, log, and shut down
        fprintf(stderr, "Failed to create SDL renderer: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        exit(EXIT_FAILURE);
    }

    // initialize SDL_image
    if ((IMG_Init(IMG_INIT_PNG) & IMG_INIT_PNG) == 0) {
        // Failed. Log an error, clean up, and exit
        fprintf(stderr, "Failed to initialize SDL_image: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        exit(EXIT_FAILURE);
    }

    // load image of blue circle (mass 1)
    SDL_Texture* image_1 = IMG_LoadTexture(renderer, "circle_blue.png");
    if (image_1 == NULL) {
        // detect error, log, and shut down
        fprintf(stderr, "Failed to create SDL texture (image): %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        exit(EXIT_FAILURE);
    }

    // load image of magenta circle (mass 2)
    SDL_Texture* image_2 = IMG_LoadTexture(renderer, "circle_magenta.png");
    if (image_2 == NULL) {
        // detect error, log, and shut down
        fprintf(stderr, "Failed to create SDL texture (image): %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        exit(EXIT_FAILURE);
    }

    int i = 0;
    double *pend1_x_arr = (double*)calloc(total_steps, sizeof(double));
    double *pend1_y_arr = (double*)calloc(total_steps, sizeof(double));
    double *pend2_x_arr = (double*)calloc(total_steps, sizeof(double));
    double *pend2_y_arr = (double*)calloc(total_steps, sizeof(double));

    // start game loop
    bool running = true;
    while (running) {
        // event handling
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                // get SDL_QUIT event when user clicks on 'x' in window
                running = false;
            }
        }

        bool running_animation = true;
        while (running_animation) {
            // end animation game loop if all pendulum positions have been printed
            if (i == total_steps) {
                running_animation = false;
                break;
            }

            // update state: clear renderer
            SDL_RenderClear(renderer);

            // put pendulum mass centers into screen coordinates array
            pend1_x_arr[i] = SCREEN_WIDTH/2 + x_1_arr[i];
            pend1_y_arr[i] = SCREEN_HEIGHT/2 + y_1_arr[i];
            pend2_x_arr[i] = SCREEN_WIDTH/2 + x_2_arr[i];
            pend2_y_arr[i] =  SCREEN_HEIGHT/2 + y_2_arr[i];

            // update display
            // fill screen with grey pixels
            SDL_SetRenderDrawColor(renderer, 200, 200, 200, 255);
            SDL_Rect screen_rect = {.x = 0, .y = 0, .w = SCREEN_WIDTH, .h = SCREEN_HEIGHT};
            SDL_RenderFillRect(renderer, &screen_rect);
       
            // draw pendulum pivot
            SDL_Rect pivot = {
                .x = SCREEN_WIDTH/2 - PIVOT_SIZE/2,
                .y = SCREEN_HEIGHT/2 - PIVOT_SIZE/2,
                .w = PIVOT_SIZE,
                .h = PIVOT_SIZE
            };
            SDL_SetRenderDrawColor(renderer, 50, 50, 50, 255);
            SDL_RenderFillRect(renderer, &pivot);

            // calculate mass image sizes based on weights
            double mass_size1 = (constants->m_1 * 4) + 10;
            double mass_size2 = (constants->m_2 * 4) + 10;

            // set pendulum mass coordinates for images
            SDL_Rect pend1 = {
                .x = pend1_x_arr[i] - mass_size1/2,
                .y = pend1_y_arr[i] - mass_size1/2,
                .w = mass_size1,
                .h = mass_size1
            };

            SDL_Rect pend2 = {
                .x = pend2_x_arr[i] - mass_size2/2,
                .y = pend2_y_arr[i] - mass_size2/2,
                .w = mass_size2,
                .h = mass_size2
            };

            // draw lines connecting the masses
            SDL_RenderDrawLine(renderer, SCREEN_WIDTH/2, SCREEN_HEIGHT/2, pend1_x_arr[i], pend1_y_arr[i]);
            SDL_RenderDrawLine(renderer, pend1_x_arr[i], pend1_y_arr[i], pend2_x_arr[i], pend2_y_arr[i]);
        
            // print pendulum masses
            SDL_RenderCopy(renderer, image_1, NULL, &pend1);
            SDL_RenderCopy(renderer, image_2, NULL, &pend2);

            // update display
            SDL_RenderPresent(renderer);
            i++;

            // pause before drawing the next frame
            SDL_Delay(1000 / FRAME_RATE);
            }

        // clear renderer
        SDL_RenderClear(renderer);

        // fill screen with grey pixels
        SDL_SetRenderDrawColor(renderer, 200, 200, 200, 255);
        SDL_Rect screen_rect = {.x = 0, .y = 0, .w = SCREEN_WIDTH, .h = SCREEN_HEIGHT};
        SDL_RenderFillRect(renderer, &screen_rect);

        // print pendulum path (mass 1)
        SDL_SetRenderDrawColor(renderer, 12, 80, 140, 255);
        for (int i = 0; i < total_steps; i++) {
            SDL_RenderDrawPoint(renderer, pend1_x_arr[i], pend1_y_arr[i]);
        }

        // print pendulum path (mass 2)
        SDL_SetRenderDrawColor(renderer, 140, 12, 80, 255);
        for (int i = 0; i < total_steps; i++) {
            SDL_RenderDrawPoint(renderer, pend2_x_arr[i], pend2_y_arr[i]);
        }

        // update display
        SDL_RenderPresent(renderer);
    }

    // shut down SDL and exit
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    // free memory
    free(init_conds);
    free(constants);
    free(state);
    free(cart);
    free(pend1_x_arr);
    free(pend1_y_arr);
    free(pend2_x_arr);
    free(pend2_y_arr);
}