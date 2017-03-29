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
    std::string line;
    std::string nameArray[500][2];
    std::string columnHeading[117];
    float values[500][115];
    std::map<std::string, int> columnAbbrvMap;
    float* selectedColumn = NULL;
    float* selectedColumnPartition = NULL;
    float** selectedColumns = NULL;

    MPI_Init(NULL, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    if(500%worldSize != 0)
    {
        if(rank == 0)
        {
            std::cout << "Error, number of process does not evenly divide 500" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if(rank == 0)
    {
        std::ifstream file( "500_Cities__City-level_Data__GIS_Friendly_Format_.csv" );
        if( file.is_open() )
        {
            getline( file, line );
            std::stringstream   columnStream( line );
            std::string         columnValue;
            char startLetter = 'A' - 1;
            int abbrvIndex = 0;
            std::string columnAbbrv = "A";

            for( int i = 0; i < 117; i++ )
            {
                columnAbbrvMap[columnAbbrv] = i;
                columnAbbrv.at( abbrvIndex ) += 1;

                if( columnAbbrv.at( abbrvIndex ) > 90 )
                {
                    if( columnAbbrv.length() == 1 )
                    {
                        columnAbbrv.at( abbrvIndex ) = 65;
                        abbrvIndex++;
                        startLetter++;
                        columnAbbrv.insert( 0, 1, startLetter );
                    }
                    else
                    {
                        startLetter++;
                        columnAbbrv.at( abbrvIndex ) = 65;
                        columnAbbrv.at( 0 ) += 1;
                    }
                }
            }

            int columnIndex = 0;
            while( getline( columnStream, columnValue, ',' ) )
            {
                columnHeading[columnIndex] = columnValue;
                columnIndex++;
            }

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
                        if( value.at( idx ) == '\"' )
                        {
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

        // selectedColumn = (float*)malloc(sizeof(float) * 500);
        // for(int i = 0; i<500; i++)
        // {
        //     selectedColumn[i] = values[i][columnAbbrvMap[argv[3]]-2];
        // }
    }

    std::string strategy = argv[1];

    if(strategy.compare("sr") == 0)
    {
        if(rank == 0)
        {
            selectedColumn = (float*)malloc(sizeof(float) * 500);
            for(int i = 0; i<500; i++)
            {
                selectedColumn[i] = values[i][columnAbbrvMap[argv[3]]-2];
            }
        }

        int elementsPerProcess = 500/worldSize;
        selectedColumnPartition = (float*)malloc(sizeof(float) * elementsPerProcess);

        MPI_Scatter(selectedColumn, elementsPerProcess, MPI_FLOAT, selectedColumnPartition,
            elementsPerProcess, MPI_FLOAT, 0, MPI_COMM_WORLD);

        std::string operation = argv[2];

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
                std::cout << "Average " << columnHeading[columnAbbrvMap[argv[3]]] << " = " << overallSum / (worldSize * elementsPerProcess) << std::endl;
            }
        }
        else if(operation.compare("min") == 0)
        {
            float processMin = selectedColumnPartition[0];

            for(int i = 1; i < elementsPerProcess; i++)
            {
                if(selectedColumnPartition[i] < processMin)
                {
                    processMin = selectedColumnPartition[i];
                }
            }

            float overallMin;

            MPI_Reduce(&processMin, &overallMin, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);

            if(rank == 0)
            {
                int minIndex = 0;
                for(minIndex = 0; minIndex < 500; minIndex++)
                {
                    if(selectedColumn[minIndex] == overallMin)
                    {
                        break;
                    }
                }

                std::cout << nameArray[minIndex][1] << ", " << nameArray[minIndex][2] << ", " << columnHeading[columnAbbrvMap[argv[3]]] << " = " << overallMin << std::endl;
            }
        }
        else if(operation.compare("max") == 0)
        {
            float processMax = selectedColumnPartition[0];

            for(int i = 1; i < elementsPerProcess; i++)
            {
                if(selectedColumnPartition[i] > processMax)
                {
                    processMax = selectedColumnPartition[i];
                }
            }

            float overallMax;

            MPI_Reduce(&processMax, &overallMax, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

            if(rank == 0)
            {
                int maxIndex = 0;
                for(maxIndex = 0; maxIndex < 500; maxIndex++)
                {
                    if(selectedColumn[maxIndex] == overallMax)
                    {
                        break;
                    }
                }

                std::cout << std::setprecision(15)
                          << nameArray[maxIndex][1] << ", " << nameArray[maxIndex][2] << ", " << columnHeading[columnAbbrvMap[argv[3]]] << " = " << overallMax << std::endl;
            }
        }
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
        if((argc - 3) != worldSize)
        {
            if(rank == 0)
            {
                std::cout << "Error, number of process does not match number of columns." << std::endl;
            }
        }
        else
        {
            selectedColumns = alloc_2d_float(worldSize, 500);

            if(rank == 0)
            {
                for(int i = 0; i < worldSize; i++)
                {
                    for(int j = 0; j < 500; j++)
                    {
                        selectedColumns[i][j] = values[j][columnAbbrvMap[argv[i+3]]-2];
                        //std::cout << selectedColumns[i][j] << std::endl;
                    }
                }
            }

            MPI_Bcast(&(selectedColumns[0][0]), worldSize * 500, MPI_FLOAT, 0, MPI_COMM_WORLD);

            std::string operation = argv[2];

            if(operation.compare("max") == 0)
            {
                float maxValue;
                float* maxValues;

                if(rank == 0)
                {
                    maxValues = (float*)malloc(sizeof(float) * worldSize);
                }

                for(int i = 0; i < worldSize; i++)
                {
                    if(rank == i)
                    {
                        maxValue = selectedColumns[i][0];

                        for(int j = 1; j < 500; j++)
                        {
                            if(selectedColumns[i][j] > maxValue)
                            {
                                maxValue = selectedColumns[i][j];
                            }
                        }
                    }
                }

                MPI_Gather(&maxValue, 1, MPI_FLOAT, maxValues, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

                if(rank == 0)
                {
                    for(int i = 0; i < worldSize; i++)
                    {
                        int index = 0;
                        for(index = 0; index < 500; index++)
                        {
                            if(selectedColumns[i][index] == maxValues[i])
                            {
                                break;
                            }
                        }

                        std::cout << operation <<  " " << columnHeading[columnAbbrvMap[argv[i+3]]] << " = " << maxValues[i]
                        << " " <<  nameArray[index][1] << " " << nameArray[index][0] << std::endl;
                    }

                    free(maxValues);
                }
            }
            else if(operation.compare("min") == 0)
            {
                float minValue;
                float* minValues;

                if(rank == 0)
                {
                    minValues = (float*)malloc(sizeof(float) * worldSize);
                }

                for(int i = 0; i < worldSize; i++)
                {
                    if(rank == i)
                    {
                        minValue = selectedColumns[i][0];

                        for(int j = 1; j < 500; j++)
                        {
                            if(selectedColumns[i][j] < minValue)
                            {
                                minValue = selectedColumns[i][j];
                            }
                        }
                    }
                }

                MPI_Gather(&minValue, 1, MPI_FLOAT, minValues, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

                if(rank == 0)
                {
                    for(int i = 0; i < worldSize; i++)
                    {
                        int index = 0;
                        for(index = 0; index < 500; index++)
                        {
                            if(selectedColumns[i][index] == minValues[i])
                            {
                                break;
                            }
                        }

                        std::cout << operation <<  " " << columnHeading[columnAbbrvMap[argv[i+3]]] << " = " << minValues[i]
                        << " " <<  nameArray[index][1] << " " << nameArray[index][0] << std::endl;
                    }

                    free(minValues);
                }
            }
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
