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

import (
	"bytes"
	"encoding/binary"
	"encoding/hex"
	"fmt"
	"unsafe"
)

type Record struct {
	*bamRecord
	filled     bool
	cigar      []CigarOp
	nameStr    string
	seqBytes   []byte
	qualScores []int8
	auxBytes   []byte
	auxTags    []Aux
}

func (self *Record) ReferenceID() int {
	self.fillData()
	return int(self.tid())
}

func (self *Record) Name() string {
	self.fillData()
	return self.nameStr
}

func (self *Record) Seq() []byte {
	self.fillData()
	return self.seqBytes
}

func (self *Record) Quality() []int8 {
	self.fillData()
	return self.qualScores
}

func (self *Record) Cigar() []CigarOp {
	self.fillData()
	return self.cigar
}

func (self *Record) Tag(tag []byte) (v Aux, ok bool) {
	self.fillData()
	for i := range self.auxTags {
		if bytes.Compare(self.auxTags[i][:2], tag) == 0 {
			return self.auxTags[i], true
		}
	}
	return
}

func (self *Record) Tags() []Aux {
	self.fillData()
	return self.auxTags
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

func (self *Record) Score() byte {
	return self.qual()
}

func (self *Record) Flags() Flags {
	return self.flag()
}

func (self *Record) Strand() int8 {
	if self.Flags()&Reverse == Reverse {
		return -1
	}
	return 1
}

var (
	bamNT16TableRev = [16]byte{'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'}
	bamNT16Table    = [256]byte{
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0x1, 0x2, 0x4, 0x8, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0x0, 0xf, 0xf,
		0xf, 0x1, 0xe, 0x2, 0xd, 0xf, 0xf, 0x4, 0xb, 0xf, 0xf, 0xc, 0xf, 0x3, 0xf, 0xf,
		0xf, 0xf, 0x5, 0x6, 0x8, 0xf, 0x7, 0x9, 0xf, 0xa, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0x1, 0xe, 0x2, 0xd, 0xf, 0xf, 0x4, 0xb, 0xf, 0xf, 0xc, 0xf, 0x3, 0xf, 0xf,
		0xf, 0xf, 0x5, 0x6, 0x8, 0xf, 0x7, 0x9, 0xf, 0xa, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
	}
)

func (self *Record) fillData() {
	if self.filled || self.bamRecord.b == nil {
		return
	}

	d := self.dataUnsafe()
	var s, e int

	// Get query name.
	s, e = 0, int(self.lQname())
	self.nameStr = string(d[s : e-1])

	// Get CIGAR data.
	nCigar := self.nCigar()
	s, e = e, e+int(nCigar<<2) // CIGAR represented as C.uint32 so length is 4*n_cigar
	self.cigar = make([]CigarOp, nCigar)
	err := binary.Read(bytes.NewBuffer(d[s:e]), endian, &self.cigar)
	if err != nil {
		panic(fmt.Sprintf("boom: binary.Read failed: %v", err))
	}

	// Get sequence data.
	lQqual := int(self.lQseq())
	lQseq := (lQqual + 1) >> 1
	// Extract nucleotide nybbles.
	s, e = e, e+lQseq
	self.seqBytes = make([]byte, lQqual)
	for i, c := range d[s:e] {
		i2 := i << 1
		self.seqBytes[i2] = bamNT16TableRev[c>>4]
		if i2++; i2 == len(self.seqBytes) {
			break
		}
		self.seqBytes[i2] = bamNT16TableRev[c&0xf]
	}
	// Get quality scores.
	s, e = e, e+lQqual
	self.qualScores = make([]int8, lQqual)
	q := d[s:e]
	copy(self.qualScores, *(*[]int8)(unsafe.Pointer(&q)))

	// Get auxilliary tags.
	lAux := self.lAux()
	s, e = e, e+lAux
	self.auxBytes = make([]byte, lAux)
	copy(self.auxBytes, d[s:e])
	self.auxTags = parseAux(self.auxBytes)

	self.filled = true
}

// Compact Idiosyncratic Gapped Alignment Report
type CigarOp uint32

func (co CigarOp) Type() CigarOpType { return CigarOpType(co & 0xf) }
func (co CigarOp) Len() int          { return int(co >> 4) }
func (co CigarOp) String() string    { return fmt.Sprintf("%d%s", co.Len(), co.Type().String()) }

type CigarOpType byte

const (
	CigarMatch CigarOpType = iota
	CigarInsertion
	CigarDeletion
	CigarSkipped
	CigarSoftClipped
	CigarHardClipped
	CigarPadded
	CigarEqual
	CigarMismatch
	lastCigar
)

var cigarOps = []string{"M", "I", "D", "N", "S", "H", "P", "=", "X", "?"}

func (ct CigarOpType) String() string {
	if ct < 0 || ct > lastCigar {
		ct = lastCigar
	}
	return cigarOps[ct]
}

type Aux []byte

var (
	jumps = [256]int{
		'A': 1,
		'c': 1, 'C': 1,
		's': 2, 'S': 2,
		'i': 4, 'I': 4,
		'f': 4,
		'Z': -1,
		'H': -1,
		'B': -1,
	}
	auxTypes = [256]byte{
		'A': 'A',
		'c': 'i', 'C': 'i',
		's': 'i', 'S': 'i',
		'i': 'i', 'I': 'i',
		'f': 'f',
		'Z': 'Z',
		'H': 'H',
		'B': 'B',
	}
)

func parseAux(aux []byte) (aa []Aux) {
	for i := 0; i+2 < len(aux); {
		t := aux[i+2]
		switch j := jumps[t]; {
		case j > 0:
			j += 3
			aa = append(aa, Aux(aux[i:i+j]))
			i += j
		case j < 0:
			switch t {
			case 'Z', 'H':
				var (
					j int
					v byte
				)
				for j, v = range aux[i:] {
					if v == 0 { // C string termination
						break // Truncate terminal zero.
					}
				}
				aa = append(aa, Aux(aux[i:i+j]))
				i += j + 1
			case 'B':
				var length int32
				err := binary.Read(bytes.NewBuffer([]byte(aux[i+4:i+8])), endian, &length)
				if err != nil {
					panic(fmt.Sprintf("boom: binary.Read failed: %v", err))
				}
				j = int(length)*jumps[aux[i+3]] + int(unsafe.Sizeof(length)) + 4
				aa = append(aa, Aux(aux[i:i+j]))
				i += j
			}
		default:
			panic(fmt.Sprintf("boom: unrecognised optional field type: %q", t))
		}
	}
	return
}

func buildAux(aa []Aux) (aux []byte) {
	for _, a := range aa {
		// TODO: validate each a
		aux = append(aux, []byte(a)...)
	}
	return
}

func (self Aux) String() string {
	return fmt.Sprintf("%s:%c:%v", []byte(self[:2]), auxTypes[self.Type()], self.Value())
}

func (self Aux) Tag() string { return string(self[:2]) }
func (self Aux) Type() byte  { return self[2] }
func (self Aux) Value() (v interface{}) {
	switch t := self.Type(); t {
	case 'A':
		return self[3]
	case 'c':
		return int8(self[3])
	case 'C':
		return uint(self[3])
	case 's':
		s := int16(0)
		err := binary.Read(bytes.NewBuffer([]byte(self[4:6])), endian, &s)
		if err != nil {
			panic(fmt.Sprintf("boom: binary.Read failed: %v", err))
		}
		return s
	case 'S':
		S := uint16(0)
		err := binary.Read(bytes.NewBuffer([]byte(self[4:6])), endian, &S)
		if err != nil {
			panic(fmt.Sprintf("boom: binary.Read failed: %v", err))
		}
		return S
	case 'i':
		i := int32(0)
		err := binary.Read(bytes.NewBuffer([]byte(self[4:8])), endian, &i)
		if err != nil {
			panic(fmt.Sprintf("boom: binary.Read failed: %v", err))
		}
		return i
	case 'I':
		I := uint32(0)
		err := binary.Read(bytes.NewBuffer([]byte(self[4:8])), endian, &I)
		if err != nil {
			panic(fmt.Sprintf("boom: binary.Read failed: %v", err))
		}
		return I
	case 'f':
		f := float32(0)
		err := binary.Read(bytes.NewBuffer([]byte(self[4:8])), endian, &f)
		if err != nil {
			panic(fmt.Sprintf("boom: binary.Read failed: %v", err))
		}
		return f
	case 'Z': // Z and H Require that parsing stops before the terminating zero.
		return string(self[3:])
	case 'H':
		h := make([]byte, hex.DecodedLen(len(self[3:])))
		_, err := hex.Decode(h, []byte(self[3:]))
		if err != nil {
			panic(fmt.Sprintf("boom: hex decoding error: %v", err))
		}
		return h
	case 'B':
		var length int32
		err := binary.Read(bytes.NewBuffer([]byte(self[4:8])), endian, &length)
		if err != nil {
			panic(fmt.Sprintf("boom: binary.Read failed: %v", err))
		}
		switch t := self[3]; t {
		case 'c':
			c := self[4:]
			return *(*[]int8)(unsafe.Pointer(&c))
		case 'C':
			return []uint8(self[4:])
		case 's':
			Bs := make([]int16, length)
			err := binary.Read(bytes.NewBuffer([]byte(self[8:])), endian, &Bs)
			if err != nil {
				panic(fmt.Sprintf("boom: binary.Read failed: %v", err))
			}
			return Bs
		case 'S':
			BS := make([]uint16, length)
			err := binary.Read(bytes.NewBuffer([]byte(self[8:])), endian, &BS)
			if err != nil {
				panic(fmt.Sprintf("boom: binary.Read failed: %v", err))
			}
			return BS
		case 'i':
			Bi := make([]int32, length)
			err := binary.Read(bytes.NewBuffer([]byte(self[8:])), endian, &Bi)
			if err != nil {
				panic(fmt.Sprintf("boom: binary.Read failed: %v", err))
			}
			return Bi
		case 'I':
			BI := make([]uint32, length)
			err := binary.Read(bytes.NewBuffer([]byte(self[8:])), endian, &BI)
			if err != nil {
				panic(fmt.Sprintf("boom: binary.Read failed: %v", err))
			}
			return BI
		case 'f':
			Bf := make([]float32, length)
			err := binary.Read(bytes.NewBuffer([]byte(self[8:])), endian, &Bf)
			if err != nil {
				panic(fmt.Sprintf("boom: binary.Read failed: %v", err))
			}
			return Bf
		default:
			panic(fmt.Sprintf("boom: unknown array type %q", t))
		}
	default:
		panic(fmt.Sprintf("boom: unknown type %q", t))
	}
	return
}
