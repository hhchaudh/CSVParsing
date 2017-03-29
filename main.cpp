// Haaris Chaudhry
// Project 2
// March 29, 2017

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <iomanip>
#include <mpi.h>

// Function obtained from http://stackoverflow.com/questions/5901476/sending-and-receiving-2d-array-over-mpi
// Returns a 2d array with contiguous memory so that the array can be successfully passed to MPI functions
float **alloc_2d_float(int rows, int cols) {
    float *data = (float *)malloc(rows*cols*sizeof(float));
    float **array= (float **)malloc(rows*sizeof(float*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

int main(int argc, char** argv)
{
    std::string                         line;   // Stores each line of the file during the file reading process
    std::string                         nameArray[500][2];  // Stores the first two columns
    std::string                         columnHeading[117]; // Stores all of the column headings
    float                               values[500][115];   // Stores the 115 columns after the columns stored by nameArray
    std::map<std::string, int>          columnAbbrvMap;     // Maps column alphabet abbreviations to specific indexes

    // The following 3 arrays will be used for MPI operations
    float*                              selectedColumn = NULL;
    float*                              selectedColumnPartition = NULL;
    float**                             selectedColumns = NULL;

    // Taken from http://mpi-forum.org/docs/mpi-1.1/mpi-11-html/node79.html
    // This struct is used in conjunction with MPI_FLOAT_INT in order to use MPI_MAXLOC AND MPI_MINLOC effectively
    struct valueInfo {
        float value;
        int index;
    } in, out;

    MPI_Init(NULL, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    // My awful parsing methodology
    if(rank == 0)
    {
        std::ifstream file( "500_Cities__City-level_Data__GIS_Friendly_Format_.csv" );
        if( file.is_open() )
        {
            getline( file, line );
            std::stringstream   columnStream( line );
            std::string         columnValue;
            char startLetter = 'A' - 1; // The start letter is used once we get to over 26 columns, it is prepended to columnAbbrv
            int abbrvIndex = 0; // The index of the letter that changes with each iteration
            std::string columnAbbrv = "A";  // The current value of the column abbreviation

            for( int i = 0; i < 117; i++ )
            {
                columnAbbrvMap[columnAbbrv] = i;    // Starts with 'A' mapping to 0, then 'B' mapping to 1....'AA' mapping to 26...
                columnAbbrv.at( abbrvIndex ) += 1;  // Increment the letter that needs to be incremented

                if( columnAbbrv.at( abbrvIndex ) > 90 ) // We went past 'Z', need to increment the startLetter and reset the letter at abbrvIndex
                {
                    if( columnAbbrv.length() == 1 )
                    {
                        columnAbbrv.at( abbrvIndex ) = 65;
                        abbrvIndex++;
                        startLetter++;
                        columnAbbrv.insert( 0, 1, startLetter );
                    }
                    else    // Not worried about triple letter columns, this parsing only handles single and double letter columns
                    {
                        startLetter++;
                        columnAbbrv.at( abbrvIndex ) = 65;
                        columnAbbrv.at( 0 ) += 1;
                    }
                }
            }

            // Read the first line, which contains the columns names
            int columnIndex = 0;
            while( getline( columnStream, columnValue, ',' ) )
            {
                columnHeading[columnIndex] = columnValue;
                columnIndex++;
            }

            // Read the rest of the file by taking each line, then storing that line
            // in a stringstream and parsing the individual line.  We throw away
            // any values that have double quotes and store 0 in the values matrix
            // for the double quoted entry.
            unsigned int i = 0;
            while( getline( file, line ) )
            {
                std::stringstream   linestream( line );
                std::string         value;

                unsigned int j = 0;
                while( getline( linestream, value, ',' ) )
                {
                    bool printVal = true;

                    for( unsigned int idx = 0; idx < value.length(); idx++ )
                    {
                        // Ignore double quoted "range" values
                        if( value.at( idx ) == '\"' )
                        {
                            // Need an offset of -2 because the first 2 columns were used by nameArray
                            // We don't have to worry about double quotes in the first two columns
                            values[i][j - 2] = 0;
                            j++;
                            printVal = false;
                            getline( linestream, value, ',' );
                            break;
                        }
                    }

                    if( printVal )
                    {
                        if( j < 2 )
                        {
                            nameArray[i][j] = value;
                        }
                        else
                        {
                            values[i][j - 2] = std::stof( value );
                        }
                        j++;
                    }
                }
                i++;
            }
            file.close();
        }
        else
        {
            std::cout << "Unable to open file";
        }
    }

    // First three cmd line entries should always contain at least the strategy, operation, and at least one column abbreviation
    // There is no error checking for incorrect column entries i.e. Entering "12" or "AAA" for argv[3]
    // Checking the min, max, or average of the first two columns will probably result in an error.

    // Here we read the first value entered in the command line, which is the strategy to use
    std::string strategy = argv[1];

    // This giant if-block handles the scatter-reduce methodology
    if(strategy.compare("sr") == 0)
    {
        // Abort if we get a number of processe that does not evenly divide 500
        if(500%worldSize != 0)
        {
            if(rank == 0)
            {
                std::cout << "Error, number of processes does not evenly divide 500" << std::endl;
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        // Rank 0 will create an array that will store all of the values, in order, of the requested column
        if(rank == 0)
        {
            selectedColumn = (float*)malloc(sizeof(float) * 500);
            for(int i = 0; i<500; i++)
            {
                // We have to use an offset for the values 2d array because the nameArray took the first two columns
                selectedColumn[i] = values[i][columnAbbrvMap[argv[3]]-2];
            }
        }

        // Each rank will create arrays that are identical in size to hold a piece of the data contained in the array created in rank 0
        int elementsPerProcess = 500/worldSize;
        selectedColumnPartition = (float*)malloc(sizeof(float) * elementsPerProcess);

        // Rank 0 scatters its selectedColumn array to be evenly distributed to all ranks, the distribution will be contained in selectedColumnPartition
        MPI_Scatter(selectedColumn, elementsPerProcess, MPI_FLOAT, selectedColumnPartition,
            elementsPerProcess, MPI_FLOAT, 0, MPI_COMM_WORLD);

        // Check for the operation that is requested (min, max, avg, or number)
        std::string operation = argv[2];

        // For average, each rank computes the sum of the values in its array and then those values are reduced to one sum using MPI_SUM
        // The overall sum is divided by the total number of elements distributed (shold be 500 every time).
        if(operation.compare("avg") == 0)
        {
            float processSum = 0;
            for(int i = 0; i < elementsPerProcess; i++)
            {
                processSum += selectedColumnPartition[i];
            }

            float overallSum;

            MPI_Reduce(&processSum, &overallSum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

            if(rank == 0)
            {
                // We get the correct columnHeading by looking at the argv[3] and getting the index associated with the column Abbreviation
                // This index corresponds to the index that the appropriate columnHeading would be found at.
                std::cout << "Average " << columnHeading[columnAbbrvMap[argv[3]]] << " = " << overallSum / (worldSize * elementsPerProcess) << std::endl;
            }
        }
        // For the min operation, we use the valueInfo struct to store the minimum value for each rank, and then MPI_MINLOC is used to reduce these structs to the one
        // that contains the smallest value.  The out struct is used by rank 0 to display the results of the reduction.  A similar process is used for finding the max
        // value.  Using a struct allows us to have the index value (which is computed by each rank).
        else if(operation.compare("min") == 0)
        {
            in.value = selectedColumnPartition[0];
            in.index = 0;

            for(int i = 1; i < elementsPerProcess; i++)
            {
                if(selectedColumnPartition[i] < in.value)
                {
                    in.value = selectedColumnPartition[i];  // Save the smallest value
                    in.index = i;                           // Save the index of this smallest value
                }
            }
            // Calculate the index of the local min with respect to the global data
            in.index = rank * elementsPerProcess + in.index;

            MPI_Reduce(&in, &out, 1, MPI_FLOAT_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);

            if(rank == 0)
            {
                float overallMin = out.value;
                int minIndex = out.index;

                std::cout << nameArray[minIndex][1] << ", " << nameArray[minIndex][2] << ", " << columnHeading[columnAbbrvMap[argv[3]]] << " = " << overallMin << std::endl;
            }
        }
        else if(operation.compare("max") == 0)
        {
            in.value = selectedColumnPartition[0];
            in.index = 0;

            for(int i = 1; i < elementsPerProcess; i++)
            {
                if(selectedColumnPartition[i] > in.value)
                {
                    in.value = selectedColumnPartition[i];
                    in.index = i;
                }
            }
            in.index = rank * elementsPerProcess + in.index;

            MPI_Reduce(&in, &out, 1, MPI_FLOAT_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);

            if(rank == 0)
            {
                float overallMax = out.value;
                int maxIndex = out.index;

                std::cout << std::setprecision(15) <<
                            nameArray[maxIndex][1] << ", " << nameArray[maxIndex][2] << ", " << columnHeading[columnAbbrvMap[argv[3]]] << " = " << overallMax << std::endl;
            }
        }
        // The number operations work similarly to average, where we get the sum from each rank and then compute the overall sum using MPI_SUM
        else if(operation.compare("number") == 0)
        {
            std::string numberOperation = argv[4];
            float comparedValue = std::stof(argv[5]);

            if(numberOperation.compare("gt") == 0)
            {
                int processValsGreater = 0;

                for(int i = 0; i < elementsPerProcess; i++)
                {
                    if(selectedColumnPartition[i] > comparedValue)
                    {
                        processValsGreater++;
                    }
                }

                int overallValsGreater;

                MPI_Reduce(&processValsGreater, &overallValsGreater, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

                if(rank == 0)
                {
                    std::cout << "Number cities with " << columnHeading[columnAbbrvMap[argv[3]]] << " gt " << comparedValue << " = " << overallValsGreater << std::endl;
                }
            }
            else if(numberOperation.compare("lt") == 0)
            {
                int processValsLess = 0;

                for(int i = 0; i < elementsPerProcess; i++)
                {
                    if(selectedColumnPartition[i] < comparedValue)
                    {
                        processValsLess++;
                    }
                }

                int overallValsLess;

                MPI_Reduce(&processValsLess, &overallValsLess, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

                if(rank == 0)
                {
                    std::cout << "Number cities with " << columnHeading[columnAbbrvMap[argv[3]]] << " lt " << comparedValue << " = " << overallValsLess << std::endl;
                }
            }
        }

        if(rank == 0)
        {
            free(selectedColumn);
        }
        free(selectedColumnPartition);
    }
    else if(strategy.compare("bg") == 0)
    {
        // If the number of columns doesn't match the number of processes to be used, then abort.
        if((argc - 3) != worldSize)
        {
            if(rank == 0)
            {
                std::cout << "Error, number of processes does not match number of columns." << std::endl;
            }
        }
        else
        {
            // Create a specialized 2d array
            selectedColumns = alloc_2d_float(worldSize, 500);

            // rank 0 populates the array
            if(rank == 0)
            {
                for(int i = 0; i < worldSize; i++)
                {
                    for(int j = 0; j < 500; j++)
                    {
                        // The arrays in selectedColumns will correspond to the columns given in the command line arguments
                        // i.e selectedColumns[0] corresponds to the first requested column
                        // selectedColumns[1] corresponds to the second requested columns, etc
                        selectedColumns[i][j] = values[j][columnAbbrvMap[argv[i+3]]-2];
                    }
                }
            }

            // Broadcast the array to all ranks
            MPI_Bcast(&(selectedColumns[0][0]), worldSize * 500, MPI_FLOAT, 0, MPI_COMM_WORLD);

            std::string operation = argv[2];

            // The max operation is similar to the max operation for scatter reduce in that we'll use the valueInfo struct to obtain
            // max values from each rank.  The difference is that each rank will be in charge of a column and will find the max
            // of the entire column rather than finding the max of a piece of one column.  Rank 0 will have an array of valueInfo structs
            // that will be used to store the max values and the index (row) that the max value was found on.  We use this index value
            // to show the corresponding city and state in the nameArray. The min operation works similarly.
            if(operation.compare("max") == 0)
            {
                struct valueInfo* maxValues = NULL;

                if(rank == 0)
                {
                    // Create an array of valueInfo structs, this will be used to store valueInfo structs when we use MPI_Gather
                    maxValues = (valueInfo*)malloc(worldSize * sizeof(*maxValues));
                }

                for(int i = 0; i < worldSize; i++)
                {
                    // Each rank value will correspond to a specific row value in the 2d array
                    // This way, when we return values using MPI_Gather, we'll maintain the correct order
                    // and can therefore derive the appropriate column headings to display.
                    if(rank == i)
                    {
                        in.value = selectedColumns[i][0];

                        for(int j = 1; j < 500; j++)
                        {
                            if(selectedColumns[i][j] > in.value)
                            {
                                in.value = selectedColumns[i][j];       // Store the max value
                                in.index = j;                           // Store index of max value
                            }
                        }
                        break;  // The rank has done its work, no need to iterate further
                    }
                }

                // Use MPI_Gather to store results of each rank in rank 0's maxValues array
                MPI_Gather(&in, 1, MPI_FLOAT_INT, maxValues, 1, MPI_FLOAT_INT, 0, MPI_COMM_WORLD);

                // rank 0 iterates through the results
                if(rank == 0)
                {
                    for(int i = 0; i < worldSize; i++)
                    {
                        // Get the index so that we can find the city that corresponds to the given max data
                        int index = maxValues[i].index;

                        std::cout << operation <<  " " << columnHeading[columnAbbrvMap[argv[i+3]]] << " = " << maxValues[i].value
                        << " " <<  nameArray[index][1] << " " << nameArray[index][0] << std::endl;
                    }

                    free(maxValues);
                }
            }
            else if(operation.compare("min") == 0)
            {
                struct valueInfo* minValues = NULL;

                if(rank == 0)
                {
                    minValues = (valueInfo*)malloc(worldSize * sizeof(*minValues));
                }

                for(int i = 0; i < worldSize; i++)
                {
                    if(rank == i)
                    {
                        in.value = selectedColumns[i][0];

                        for(int j = 1; j < 500; j++)
                        {
                            if(selectedColumns[i][j] < in.value)
                            {
                                in.value = selectedColumns[i][j];
                                in.index = j;
                            }
                        }
                        break;
                    }
                }

                MPI_Gather(&in, 1, MPI_FLOAT_INT, minValues, 1, MPI_FLOAT_INT, 0, MPI_COMM_WORLD);

                if(rank == 0)
                {
                    for(int i = 0; i < worldSize; i++)
                    {
                        int index = minValues[i].index;

                        std::cout << operation <<  " " << columnHeading[columnAbbrvMap[argv[i+3]]] << " = " << minValues[i].value
                        << " " <<  nameArray[index][1] << " " << nameArray[index][0] << std::endl;
                    }

                    free(minValues);
                }
            }
            // Similar to max/min bg strategies, except that we don't need to keep track of an index and therefore can use MPI_FLOAT rather than MPI_FLOAT_INT
            else if(operation.compare("avg") == 0)
            {
                float avgValue = 0;
                float* avgValues;

                if(rank == 0)
                {
                    avgValues = (float*)malloc(sizeof(float) * worldSize);
                }

                for(int i = 0; i < worldSize; i++)
                {
                    if(rank == i)
                    {
                        float sum = 0;

                        for(int j = 0; j < 500; j++)
                        {
                            sum += selectedColumns[i][j];
                        }

                        avgValue = sum/500.0;
                        break;
                    }
                }

                MPI_Gather(&avgValue, 1, MPI_FLOAT, avgValues, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

                if(rank == 0)
                {
                    for(int i = 0; i < worldSize; i++)
                    {
                        std::cout << operation << " " << columnHeading[columnAbbrvMap[argv[i+3]]] << " = " << avgValues[i] << std::endl;
                    }

                    free(avgValues);
                }
            }

            free(selectedColumns[0]);
            free(selectedColumns);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
