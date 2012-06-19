// Copyright ©2012 Dan Kortschak <dan.kortschak@adelaide.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

package boom

// BuildIndex builds a BAM index file, filename.bai, from a sorted BAM file, filename.
// It returns any error that occured.
func BuildIndex(file string) (err error) {
	_, err = bamIndexBuild(file)
	return
}

// An Index represents an in memory BAM index.
type Index struct {
	*bamIndex
}

// LoadIndex loads a BAM index file, and returns the index in i if no error occurred.
// If an error occurred i is returned nil and the error is returned.
func LoadIndex(file string) (i *Index, err error) {
	bi, err := bamIndexLoad(file)
	if err == nil {
		i = &Index{bi}
	}
	return
}
