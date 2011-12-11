#ifndef LOCATION_H_
#define LOCATION_H_

#include <ostream>

/*
    struct for representing locations in the grid.
*/
struct Location
{
	int row, col;

	Location()
	{
		row = col = 0;
	}

	Location(int r, int c)
	{
		row = r;
		col = c;
	}

	bool operator==(const Location & o) const {
		return row == o.row && col == o.col;
	}
	bool operator!=(const Location & o) const {return !(*this == o);}
};

inline std::ostream & operator<<(std::ostream & out, const Location & loc) {
	return out << loc.row << ", " << loc.col;
}

#endif //LOCATION_H_
