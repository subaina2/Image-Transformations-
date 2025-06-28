#include <iostream>
#include<fstream>
#include <sstream>
#include<vector>
#include<cstring>
#include<string>
#include <cmath>
#include<algorithm>
#include <cstdlib>
using namespace std;

struct Image {
    char ImageFileName[100] = { '\0' };
    vector<vector<int>> ImageData;
    int cols = 400, rows = 400, maxGray = 256;
    vector<char> comment;

    bool imageLoaded;
    bool imageModified;

   
    int loadImage(const char ImageName[]) {

        ifstream FCIN("C:\\Users\\my computer & laptop\\Downloads\\khan.pgm");
        ifstream secondimage("C:\\Users\\my computer & laptop\\Downloads\\imranwithjem.pgm");
        if (!FCIN.is_open())
            return -1;

        char MagicNumber[5];
        char Comment[100];

        FCIN.getline(MagicNumber, 4);
        FCIN.getline(Comment, 100);
        FCIN >> cols >> rows >> maxGray;

        ImageData.clear();
        ImageData.resize(rows, vector<int>(cols, 0));

        for (int r = 0; r < rows; r++)
            for (int c = 0; c < cols; c++)
                FCIN >> ImageData[r][c];

        if (FCIN.fail())
            return -2;

        FCIN.close();
        imageLoaded = true;
        imageModified = false;
        strcpy_s(ImageFileName, ImageName);
        return 0;
    }

    int saveImage(const char ImageName[]) {
        ofstream FCOUT("C:\\Users\\my computer & laptop\\Downloads\\khan.pgm");
        if (!FCOUT.is_open())
            return -1;

        FCOUT << "P2\n# This is a comment\n"
            << cols << " " << rows << endl << maxGray << endl;
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++)
                FCOUT << ImageData[r][c] << " ";
            FCOUT << endl;
        }
        FCOUT.close();
        imageModified = false;
        return 0;
    }
    void verticalFlip() {
        for (int r = 0; r < rows / 2; r++)
            for (int c = 0; c < cols; c++) {
                int Temp = ImageData[r][c];
                ImageData[r][c] = ImageData[rows - r - 1][c];
                ImageData[rows - r - 1][c] = Temp;
            }
        return;
    }
    void horizontalFlip() {
        for (int r = 0; r < rows; r++)
            for (int c = 0; c < cols / 2; c++) {
                int Temp = ImageData[r][c];
                ImageData[r][c] = ImageData[r][cols - c - 1];
                ImageData[r][cols - c - 1] = Temp;
            }
        return;
    }

    void Brightness(double factor) {
        for (int r = 0; r < rows; r++)
            for (int c = 0; c < cols; c++) {
                ImageData[r][c] = static_cast<int>(ImageData[r][c] * factor);
                if (ImageData[r][c] > maxGray)
                    ImageData[r][c] = maxGray;
            }
    }

    void scaleImage(int scaleFactor) {
        int newRows = round(rows * scaleFactor);
        int newCols = round(cols * scaleFactor);

        // Allocate memory for newImageData on the heap
        vector<vector<int>> newImageData(newRows, vector<int>(newCols, 0));

        for (int r = 0; r < newRows; r++) {
            for (int c = 0; c < newCols; c++) {
                int originalR = static_cast<int>(round(r / scaleFactor));
                int originalC = static_cast<int>(round(c / scaleFactor));

                // Ensure the indices are within bounds
                originalR = max(0, min(originalR, rows - 1));
                originalC = max(0, min(originalC, cols - 1));

                newImageData[r][c] = ImageData[originalR][originalC];
            }
        }

        // Update the image data with the scaled data
        ImageData = newImageData;
        rows = newRows;
        cols = newCols;

        imageModified = true;
        return;
    }

    // Function to translate the image by dx and dy
    void Translate(int dx, int dy) {
        unsigned int newimagedata[400][400];
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                newimagedata[r][c] = 0;
            }
        }
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                int newr = r - dy;
                int newc = c - dx;

                // Check if the translated position is within bounds
                if (newr >= 0 && newr < rows && newc >= 0 && newc < cols) {
                    newimagedata[r][c] = ImageData[newr][newc];
                }
            }
        }

        // Copy the translated image back to the original image
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                ImageData[r][c] = newimagedata[r][c];
            }
        }
        imageModified = true;
        return;
    }
    void Rotate(double theta) {

        char rotationDirection;
        cout << "Enter 'a' for clockwise or 'b' for anticlockwise: ";
        cin >> rotationDirection;

        bool clockwise = (rotationDirection == 'a' || rotationDirection == 'A');

        double ct, st;
        theta *= 3.1415;
        theta /= 180;

        ct = cos(theta);
        st = sin(theta);

        vector<vector<int>> DV(rows, vector<int>(cols, 255));

        int counter = 0;
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                int nr, nc;
                nr = (round)(ct * r - st * c);
                nc = (round)(st * r + ct * c);
                if ((0 <= nr) && (nr < rows) && (0 <= nc) && (nc < cols))
                    DV[nr][nc] = ImageData[r][c];
            }
        }
        ImageData = DV;
        return;
    }


    // Function to adjust the sharpness of the image
    void Sharpen(double strength) {
        // Laplacian filter for sharpening
        vector<vector<int>> laplacianFilter = {
            {-1, -1, -1},
            {-1,  8, -1},
            {-1, -1, -1}
        };

        int filterSize = laplacianFilter.size();
        int filterRadius = filterSize / 2;

        // Create a temporary copy of the image
        vector<vector<int>> tempImageData = ImageData;

        for (int r = filterRadius; r < rows - filterRadius; r++) {
            for (int c = filterRadius; c < cols - filterRadius; c++) {
                int sum = 0;

                // Apply Laplacian filter
                for (int i = -filterRadius; i <= filterRadius; i++) {
                    for (int j = -filterRadius; j <= filterRadius; j++) {
                        sum = sum + laplacianFilter[i + filterRadius][j + filterRadius] * ImageData[r + i][c + j];
                    }
                }

                // Adjust pixel value based on the strength parameter
                int newValue = static_cast<int>(ImageData[r][c] + strength * sum);

                // Ensure the new pixel value is within the valid range [0, maxGray]
                newValue = max(0, min(newValue, maxGray));

                tempImageData[r][c] = newValue;
            }
        }

        // Update the image data with the sharpened data
        ImageData = tempImageData;

        imageModified = true;
    }

    // Function to enhance the contrast of the image using linear contrast stretching
    void EnhanceContrast() {
        int minIntensity = maxGray;
        int maxIntensity = 0;

        // Find the minimum and maximum pixel values in the image
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                minIntensity = min(minIntensity, ImageData[r][c]);
                maxIntensity = max(maxIntensity, ImageData[r][c]);
            }
        }
        // Apply linear contrast stretching to scale pixel values to the full range
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                ImageData[r][c] = static_cast<int>(
                    (ImageData[r][c] - minIntensity) * (maxGray - 1) / max(1, maxIntensity - minIntensity)
                    );
            }
        }
        imageModified = true;
    }
    void ConvertToBinary() {
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                // If the pixel value is greater than or equal to half of maxGray, set it to maxGray, otherwise set it to 0
                ImageData[r][c] = (ImageData[r][c] >= maxGray / 2) ? maxGray : 0;
            }
        }

        imageModified = true;
    }
    void ResizeImage(int newRows, int newCols) {
        // Allocate memory for newImageData on the heap
        vector<vector<int>> newImageData(newRows, vector<int>(newCols, 0));

        double scaleRow = static_cast<double>(rows) / newRows;
        double scaleCol = static_cast<double>(cols) / newCols;

        for (int r = 0; r < newRows; r++) {
            for (int c = 0; c < newCols; c++) {
                int originalR = static_cast<int>(r * scaleRow);
                int originalC = static_cast<int>(c * scaleCol);

                // Ensure the indexs are within bounds
                originalR = max(0, min(originalR, rows - 1));
                originalC = max(0, min(originalC, cols - 1));

                newImageData[r][c] = ImageData[originalR][originalC];
            }
        }

        // Update the image data with the resized data
        ImageData = newImageData;
        rows = newRows;
        cols = newCols;
        imageModified = true;
    }
    void Crop(int startRow, int startCol, int endRow, int endCol) {
        
        // Validate the input coordinates
        if (startRow < 0 || startCol < 0 || endRow >= rows || endCol >= cols || startRow >= rows || startCol >= cols || endRow < startRow || endCol < startCol) {
            cout << "Invalid cropping coordinates." << endl;
            return;
        }

        // Calculate the size of the cropped region
        int newRows = endRow - startRow + 1;
        int newCols = endCol - startCol + 1;

        // Allocate memory for newImageData on the heap
        vector<vector<int>> newImageData(newRows, vector<int>(newCols, 0));

        for (int r = 0; r < newRows; r++) {
            for (int c = 0; c < newCols; c++) {
                newImageData[r][c] = ImageData[startRow + r][startCol + c];
            }
        }

        // Update the image data with the cropped data
        ImageData = newImageData;
        rows = newRows;
        cols = newCols;
        imageModified = true;
    }
    void ManipulateIntensity(double intensityFactor) {
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                ImageData[r][c] = static_cast<int>(ImageData[r][c] * intensityFactor);

                if (ImageData[r][c] > maxGray) {
                    ImageData[r][c] = maxGray;
                }
            }
        }
        imageModified = true;
    }

    // Function to compute the derivative of the image using derivative mask filters
    void FindDerivative() {
        // Derivative mask filter
        vector<vector<int>> derivativefilter = {
            {-1, 0, 1},
            {-1, 0, 1},
            {-1, 0, 1}
        };

        int filterSize = derivativefilter.size();
        int filterRadius = filterSize / 2;

        // Create a temporary copy of the image
        vector<vector<int>> tempImageData = ImageData;

        for (int r = filterRadius; r < rows - filterRadius; r++) {
            for (int c = filterRadius; c < cols - filterRadius; c++) {
                int sum = 0;

                // Apply derivative filter
                for (int i = -filterRadius; i <= filterRadius; i++) {
                    for (int j = -filterRadius; j <= filterRadius; j++) {
                        sum = sum + derivativefilter[i + filterRadius][j + filterRadius] * ImageData[r + i][c + j];
                    }
                }

                // Set the pixel value to the computed derivative
                tempImageData[r][c] = sum;
            }
        }

        // Update the image data with the computed derivative
        ImageData = tempImageData;
        imageModified = true;
    }

    // Function to find edges in the image using derivative mask filters
    void FindEdges() {
        // Derivative mask filters for horizontal and vertical edges
        vector<vector<int>> horizontalFilter = {
            {-1, 0, 1},
            {-2, 0, 2},
            {-1, 0, 1}
        };

        vector<vector<int>> verticalFilter = {
            {-1, -2, -1},
            { 0,  0,  0},
            { 1,  2,  1}
        };

        int filterSize = horizontalFilter.size();
        int filterRadius = filterSize / 2;

        // Create a temporary copy of the image
        vector<vector<int>> tempImageData = ImageData;

        for (int r = filterRadius; r < rows - filterRadius; r++) {
            for (int c = filterRadius; c < cols - filterRadius; c++) {
                int sumX = 0, sumY = 0;

                // Apply horizontal and vertical filters
                for (int i = -filterRadius; i <= filterRadius; i++) {
                    for (int j = -filterRadius; j <= filterRadius; j++) {
                        sumX += horizontalFilter[i + filterRadius][j + filterRadius] * ImageData[r + i][c + j];
                        sumY += verticalFilter[i + filterRadius][j + filterRadius] * ImageData[r + i][c + j];
                    }
                }

                // Compute the gradient magnitude
                int gradientMagnitude = static_cast<int>(sqrt(sumX * sumX + sumY * sumY));

                // Set the pixel value to the computed gradient magnitude
                tempImageData[r][c] = gradientMagnitude;
            }
        }

        // Update the image data with the computed edges
        ImageData = tempImageData;
        imageModified = true;
    }
    // Function to combine two images side by side
    void CombineImagesSideBySide(const Image& secondImage) {
        if (rows != secondImage.rows) {
            cout << "Error: Images must have the same number of rows for combining." << endl;
            return;
        }

        int newCols = cols + secondImage.cols;

        // Allocate memory for the combined image data on the heap
        vector<vector<int>> combinedImageData(rows, vector<int>(newCols, 0));

        // Copy the data from the first image to the combined image
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                combinedImageData[r][c] = ImageData[r][c];
            }
        }

        // Copy the data from the second image to the combined image
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < secondImage.cols; c++) {
                combinedImageData[r][cols + c] = secondImage.ImageData[r][c];
            }
        }
        // Update the image data with the combined data
        ImageData = combinedImageData;
        cols = newCols;
        imageModified = true;
    }

    // Function to combine two images top to down
    void CombineImagesTopToDown(const Image& secondImage) {
        if (cols != secondImage.cols) {
            cout << "Error: Images must have the same number of columns for combining." << endl;
            return;
        }
        int newRows = rows + secondImage.rows;

        // Allocate memory for the combined image data on the heap
        vector<vector<int>> combinedImageData(newRows, vector<int>(cols, 0));

        // Copy the data from the first image to the combined image
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                combinedImageData[r][c] = ImageData[r][c];
            }
        }

        // Copy the data from the second image to the combined image
        for (int r = 0; r < secondImage.rows; r++) {
            for (int c = 0; c < cols; c++) {
                combinedImageData[rows + r][c] = secondImage.ImageData[r][c];
            }
        }

        // Update the image data with the combined data
        ImageData = combinedImageData;
        rows = newRows;
        imageModified = true;
    }


    // Function to apply linear filter to an image
    void LinearFilter(const string& imagePath, const string& filterPath) {
        // Load the image
        ifstream imageFile("C:\\Users\\my computer & laptop\\Downloads\\khan.pgm");
        if (!imageFile.is_open()) {
            cout << "Error opening image file." << endl;
            return;
        }

        // Read image properties
        string format;
        int imageSize, maxPixelValue;
        imageFile >> format >> imageSize >> imageSize >> maxPixelValue;

        // Load the image data
        vector<vector<int>> imageData(imageSize, vector<int>(imageSize));
        for (int i = 0; i < imageSize; ++i) {
            for (int j = 0; j < imageSize; ++j) {
                imageFile >> imageData[i][j];
            }
        }
        imageFile.close();

        // Load the filter
        ifstream filterFile("C:\\Users\\my computer & laptop\\Downloads\\itisfilterfile\\filterfile");
        if (!filterFile.is_open()) {
            cout << "Error opening filter file." << endl;
            return;
        }

        // Read filter data
        vector<vector<int>> filter(imageSize, vector<int>(imageSize));
        for (int i = 0; i < imageSize; ++i) {
            for (int j = 0; j < imageSize; ++j) {
                filterFile >> filter[i][j];
            }
        }
        filterFile.close();

        // Apply the filter to the image
        vector<vector<int>> newImageData(imageSize, vector<int>(imageSize, 0));

        for (int i = 1; i < imageSize - 1; ++i) {
            for (int j = 1; j < imageSize - 1; ++j) {
                int sum = 0;
                for (int x = 0; x < 3; ++x) {
                    for (int y = 0; y < 3; ++y) {
                        sum += imageData[i - 1 + x][j - 1 + y] * filter[x][y];
                    }
                }
                newImageData[i][j] = sum;
            }
        }

        // Save the changes to the original image file
        ofstream resultFile("C:\\Users\\my computer & laptop\\Downloads\\khan.pgm");
        if (!resultFile.is_open()) {
            cout << " Error saving changes to image file." << endl;
            return;
        }

        // Save the new image properties
        resultFile << format << endl;
        resultFile << imageSize << " " << imageSize << endl;
        resultFile << maxPixelValue << endl;

        // Save the new image data
        for (int i = 0; i < imageSize; ++i) {
            for (int j = 0; j < imageSize; ++j) {
                resultFile << imageData[i][j] << " ";
            }
            resultFile << endl;
        }

        resultFile.close();
    }
    void MeanFilter() {
        vector<vector<int>> newImageData(rows, vector<int>(cols, 0));

        for (int r = 1; r < rows - 1; ++r) {
            for (int c = 1; c < cols - 1; ++c) {
                int sum = 0;
                for (int i = -1; i <= 1; ++i) {
                    for (int j = -1; j <= 1; ++j) {
                        sum = sum + ImageData[r + i][c + j];
                    }
                }
                newImageData[r][c] = sum / 9; // 3x3 mean filter
            }
        }

        // Update the image data with the filtered data
        ImageData = newImageData;

        imageModified = true;
    }

    void MedianFilter() {
        vector<vector<int>> newImageData(rows, vector<int>(cols, 0));

        for (int r = 1; r < rows - 1; r++) {
            for (int c = 1; c < cols - 1; c++) {
                vector<int> value;

                for (int i = -1; i <= 1; ++i) {
                    for (int j = -1; j <= 1; ++j) {
                        value.push_back(ImageData[r + i][c + j]);
                    }
                }

                // Sort the values to find the median
                sort(value.begin(), value.end());

                newImageData[r][c] = value[4]; // 3x3 median filter
            }
        }

        // Update the image data with the filtered data
        ImageData = newImageData;
        imageModified = true;
    }
};

struct Menu {
    vector<string> menuItems;

    Menu(const char menuFile[]) {
        loadMenu(menuFile);
    }

    int loadMenu(const char menuFile[]) {
        string menufile = "C:\\Users\\my computer & laptop\\Downloads\\menufile";

        ofstream OUT;
        OUT.open("C:\\Users\\my computer & laptop\\Downloads\\menu.txt");

        if (!OUT.is_open()) {
            cerr << "Error: Unable to create menu file." << endl;
            return -1;
        }

        OUT << "1: Load Image\n2: Save Image\n3: Vertical Flip\n4: Horizontal Flip\n5: Change Brightness\n6: Scale Image\n7: Translate Image\n8: Rotate Image\n9: Sharpen the Image\n10: Enhance the Image Contrast \n11: Convert to Binary\n12: Resize the Image\n13: Crop an Imgae\n14: Manipulate image Intensity\n15: Find Image Derivatives\n16: Find Edges\n17: Combine Images\n18: Apply Linear Filter\n19: Apply Mean Filter\n20: Apply Median Filter\n21: Exit" << endl;
        OUT.close();

        ifstream IN(menuFile);
        if (!IN.is_open()) {
            cerr << "Error: Unable to open menu file for reading." << endl;
            return -1;
        }

        string line;
        while (getline(IN, line)) {
            cout << line << "  " << endl;
        }

        char menuItem[100], TotalItems[10];
        int Choices;

        IN.getline(TotalItems, 8);
        Choices = atoi(TotalItems);

        for (int i = 1; i <= Choices; i++) {
            IN.getline(menuItem, 99);
            menuItems.push_back(menuItem);
        }

        IN.close();
        return Choices;
    }

    int presentMenu() {
        int userChoice;
        size_t totalChoices = menuItems.size();

        do {
            int k = 22;
            for (size_t i = 0; i < totalChoices; i++) {
                if (menuItems[i][0] != '*') {
                    cout << k << "\t" << menuItems[i] << endl;
                    k++;
                }
            }

            cout << " Enter Your Choice (1 - " << k - 1 << " ) : ";
            cin >> userChoice;
            return userChoice;


            if (cin.fail()) {
                cin.clear();
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
                cout << "Invalid input. Please enter a number." << endl;
            }
        } while (userChoice != 22);

    }
};



int main() {
    char MenuFile[] = "C:\\Users\\my computer & laptop\\Downloads\\menu.txt";
    Image images[2];
    int activeImage = 0;
    int errorCode = 0;
    int userChoice = 0;
    int totalChoices = 22;
    Menu menu("C:\\Users\\my computer & laptop\\Downloads\\menu.txt");
    totalChoices = menu.menuItems.size();
    // Paths to image files
    string firstImagePath = "C:\\Users\\my computer & laptop\\Downloads\\khan.pgm";
    string secondImagePath = "C:\\Users\\my computer & laptop\\Downloads\\imranwithjem.pgm";
    string combinedImagePath = "C:\\Users\\my computer & laptop\\Downloads\\combined\\newimage\\";


    do {
        userChoice = menu.presentMenu();
        if (1 == userChoice) {
            char ImageFileName[100];
            cout << " Specify File Name: ";
            cin >> ImageFileName;
            errorCode = images[activeImage].loadImage(ImageFileName);
            if (errorCode == 0) {
                cout << " File Loaded Successfully " << endl;
            }
            else {
                cout << "Load Error: Code " << errorCode << endl;
            }
        }
        else if (2 == userChoice) {
            char ImageFileName[100];
            cout << " Specify File Name ";
            cin >> ImageFileName;
            errorCode = images[activeImage].saveImage(ImageFileName);
            if (errorCode == 0) {
                cout << " File Saved as " << ImageFileName << endl;
            }
            else {
                cout << "Save Error: Code " << errorCode << endl;
            }
        }
        else if (3 == userChoice) {
            images[activeImage].verticalFlip();
            cout << " Image flipped Vertically.You need to save the changes " << endl;
        }
        else if (4 == userChoice) {
            images[activeImage].horizontalFlip();
            cout << " Image Flipped Horizontally.You need to save the changes " << endl;
        }
        else if (5 == userChoice) {
            double x;
            cout << " Enter brightness level ( less than 5 to see visible change in brightness ): ";
            cin >> x;
            images[activeImage].Brightness(x);
            cout << " Brightness of the image has been changed. You need to save the changes " << endl;

        }
        else if (6 == userChoice) {
            int scaleFactor;
            cout << " Enter the scale factor: ";
            cin >> scaleFactor;
            images[activeImage].scaleImage(scaleFactor);
            cout << " Image scaled. You need to save the changes.\n";
        }
        else if (7 == userChoice)
        {
            int dx, dy;
            cout << " Enter the translation in x direction: ";
            cin >> dx;
            cout << " Enter the translation in y direction: ";
            cin >> dy;
            images[activeImage].Translate(dx, dy);
            cout << " Image translated.You need to save the changes.\n";
        }
        else if (8 == userChoice) {
            double angle;
            cout << " Enter the rotation angle in degrees: ";
            cin >> angle;
            images[activeImage].Rotate(angle);
            cout << " Image rotated. You need to save the changes.\n";

        }
        else if (9 == userChoice) {
            double sharpeningStrength;
            cout << " Enter the sharpening strength (-1.0 to 1.0): ";
            cin >> sharpeningStrength;
            images[activeImage].Sharpen(sharpeningStrength);
            cout << " Image has been sharpened.You need to save the changes" << endl;
        }
        else if (10 == userChoice) {
            images[activeImage].EnhanceContrast();
            cout << " Image has been enhanced.You need to save the changes " << endl;
        }
        else if (11 == userChoice) {
            images[activeImage].ConvertToBinary();
            cout << " Image converted to binary. You need to save the changes.\n";
        }
        else if (12 == userChoice) {
            int newRows, newCols;
            cout << " Enter the new number of rows: ";
            cin >> newRows;
            cout << " Enter the new number of columns: ";
            cin >> newCols;

            images[activeImage].ResizeImage(newRows, newCols);
            cout << " Image resized. You need to save the changes.\n";
        }
        else if (13 == userChoice) {
            int startRow, startCol, endRow, endCol;
            cout << " Enter the starting row for cropping: ";
            cin >> startRow;
            cout << " Enter the starting column for cropping: ";
            cin >> startCol;
            cout << " Enter the ending row for cropping: ";
            cin >> endRow;
            cout << " Enter the ending column for cropping: ";
            cin >> endCol;
            images[activeImage].Crop(startRow, startCol, endRow, endCol);
            cout << " Image cropped. You need to save the changes.\n";
        }
        else if (14 == userChoice) {
            double intensityFactor;
            cout << "Enter the intensity factor (e.g., 1.5 for increasing intensity, 0.5 for decreasing intensity): ";
            cin >> intensityFactor;

            images[activeImage].ManipulateIntensity(intensityFactor);
            cout << "Image intensity manipulated. You need to save the changes.\n";

        }

        else if (15 == userChoice) {
            images[activeImage].FindDerivative();
            cout << " Image derivative Calculated. You need to save the changes.\n";


        }
        else if (16 == userChoice) {
            images[activeImage].FindEdges();
            cout << "Edges found in the image. You need to save the changes.\n";
        }
        else if (17 == userChoice) {

            Image secondImage;
            errorCode = secondImage.loadImage("C:\\Users\\my computer & laptop\\Downloads\\imranwithjem.pgm");
            if (errorCode == 0) {
                int combineChoice;
                cout << "How would you like to combine the images?" << endl;
                cout << "1. Side by side" << endl;
                cout << "2. Top to down" << endl;
                cout << "Enter your choice: ";
                cin >> combineChoice;

                switch (combineChoice) {
                case 1:
                    images[activeImage].CombineImagesSideBySide(secondImage);
                    cout << " Images Combined.You need to save the changes" << endl;
                    break;
                case 2:
                    images[activeImage].CombineImagesTopToDown(secondImage);
                    cout << " Images Combined.You need to save the changes" << endl;
                    break;

                default:
                    cout << " Invalid choice. Please enter 1 or 2." << endl;
                }
            }
            else {
                cout << " Error loading the second image." << endl;
            }
        }

        else if (18 == userChoice) {

            string imagePath = "C:\\Users\\my computer & laptop\\Downloads\\khan.pgm";
            string filterPath = "C:\\Users\\my computer & laptop\\Downloads\\itisfilterfile\\filterfile";

            images[activeImage].LinearFilter(imagePath, filterPath);
            cout << "Linear filter applied. You need to save the changes.\n";
        }
        else if (19 == userChoice) {
            images[activeImage].MeanFilter();
            cout << "Mean filter applied. You need to save the changes.\n";
        }
        else if (20 == userChoice) {
            images[activeImage].MedianFilter();
            cout << "Median filter applied. You need to save the changes.\n";
        }
        else if (21 == userChoice) {
            break;

        }

    } while (userChoice != totalChoices - 1);

    return 0;
}
