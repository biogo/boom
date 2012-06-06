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

type SAMFile struct {
	*samFile
}

func OpenSAM(filename, ref string) (s *SAMFile, err error) {
	h := textHeader(ref)
	sf, err := samOpen(filename, "r", h)
	if err != nil {
		return
	}
	return &SAMFile{sf}, nil
}

var tWModes = [2]string{"w", "wh"}

func CreateSAM(filename string, ref *Header) (s *SAMFile, err error) {
	var mode string
	if ref == nil {
		mode = tWModes[0]
	} else {
		mode = tWModes[1]
	}
	sf, err := samOpen(filename, mode, ref.bamHeader)
	if err != nil {
		return
	}
	return &SAMFile{sf}, nil
}

func (self *SAMFile) Close() error {
	return self.samClose()
}

func (self *SAMFile) Read() (r *Record, n int, err error) {
	n, br, err := self.samRead()
	r = &Record{br}
	return
}

func (self *SAMFile) Write(r *Record) (n int, err error) {
	return self.samWrite(r.bamRecord)
}
