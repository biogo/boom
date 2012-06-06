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

import ()

type Record struct {
	*bamRecord
}

func (self *Record) Start() int {
	return int(self.pos())
}

func (self *Record) Len() int {
	return int(self.lQseq())
}

func (self *Record) End() int {
	return int(self.pos() + self.lQseq())
}

func (self *Record) Quality() byte {
	return self.qual()
}

func (self *Record) Flags() uint {
	return self.flag()
}
