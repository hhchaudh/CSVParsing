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

    // Taken from http://mpi-forum.org/docs/mpi-1.1/mpi-11-html/node79.html
    struct valueInfo {
        float value;
        int index;
    } in, out;

    MPI_Init(NULL, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    if(500%worldSize != 0)
    {
        if(rank == 0)
        {
            std::cout << "Error, number of processes does not evenly divide 500" << std::endl;
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
            in.value = selectedColumnPartition[0];
            in.index = 0;

            for(int i = 1; i < elementsPerProcess; i++)
            {
                if(selectedColumnPartition[i] < in.value)
                {
                    in.value = selectedColumnPartition[i];
                    in.index = i;
                }
            }
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
                std::cout << "Error, number of processes does not match number of columns." << std::endl;
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
                    }
                }
            }

            MPI_Bcast(&(selectedColumns[0][0]), worldSize * 500, MPI_FLOAT, 0, MPI_COMM_WORLD);

            std::string operation = argv[2];

            if(operation.compare("max") == 0)
            {
                struct valueInfo* maxValues = NULL;

                if(rank == 0)
                {
                    maxValues = (valueInfo*)malloc(worldSize * sizeof(*maxValues));
                }

                for(int i = 0; i < worldSize; i++)
                {
                    if(rank == i)
                    {
                        in.value = selectedColumns[i][0];

                        for(int j = 1; j < 500; j++)
                        {
                            if(selectedColumns[i][j] > in.value)
                            {
                                in.value = selectedColumns[i][j];
                                in.index = j;
                            }
                        }
                        break;
                    }
                }

                MPI_Gather(&in, 1, MPI_FLOAT_INT, maxValues, 1, MPI_FLOAT_INT, 0, MPI_COMM_WORLD);

                if(rank == 0)
                {
                    for(int i = 0; i < worldSize; i++)
                    {
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
