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

type BAMFile struct {
	*samFile
}

func OpenBAM(filename string) (b *BAMFile, err error) {
	sf, err := samOpen(filename, "rb", nil)
	if err != nil {
		return
	}
	return &BAMFile{sf}, nil
}

var bWModes = [2]string{"wb", "wbu"}

func CreateBAM(filename string, ref Header, comp bool) (b *BAMFile, err error) {
	var mode string
	if comp {
		mode = bWModes[0]
	} else {
		mode = bWModes[1]
	}
	sf, err := samOpen(filename, mode, ref.bamHeader)
	if err != nil {
		return
	}
	return &BAMFile{sf}, nil
}

func (self *BAMFile) Close() error {
	return self.samClose()
}

func (self *BAMFile) Read() (r *Record, n int, err error) {
	n, br, err := self.samRead()
	r = &Record{br}
	return
}

func (self *BAMFile) Write(r *Record) (n int, err error) {
	return self.samWrite(r.bamRecord)
}
