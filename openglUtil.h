extern float movementX;
extern float movementY;
extern float movementZ;

extern float rotateAngle;
extern float azimuthAngle;

void handleKeypress(unsigned char key, //The key that was pressed
                    int x, int y //The current mouse coordinates
					);

void handleSpecialKeyPress(int key, int x, int y);

void handleResize(int w, int h);

void setColors();

void initRendering(int Width, int Height);