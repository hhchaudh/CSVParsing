#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>

int main()
{
    std::string line;
    std::string nameArray[500][2];
    std::string columnHeading[117];
    float values[500][115];
    std::map<std::string, int> columnAbbrvMap;

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
    }
    else
    {
        std::cout << "Unable to open file";
    }

    for( unsigned int z = 0; z < 117; z++ )
    {
        if( z < 2 )
        {
            std::cout << nameArray[253][z] << std::endl;
        }
        else
        {
            std::cout << values[253][z - 2] << std::endl;
        }
    }

    std::string cat = "cat";
    cat.at( 0 ) = 65;

    std::cout << cat << std::endl;

    return 0;
}
