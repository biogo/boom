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

type Record struct {
	*bamRecord
	filled   bool
	cigarStr string
	nameStr  string
	seqStr   string
	qualStr  string
	auxStr   string
}

func (self *Record) ID() string {
	self.fillData()
	return self.nameStr
}

func (self *Record) Aux() string {
	self.fillData()
	return self.auxStr
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

func (self *Record) Strand() int8 {
	if self.Flags()&Reverse == Reverse {
		return -1
	}
	return 1
}

func (self *Record) fillData() (n int) {
	if self.filled {
		return
	}

	d := self.dataUnsafe()
	var s, e int

	s, e = 0, int(self.lQname())
	self.nameStr = string(d[s : e-1])
	self.auxStr = string(d[len(d)-int(self.lAux()):])
	// TODO: Not clear how CIGAR, seq and qual are encoded.

	self.filled = true
	return
}
