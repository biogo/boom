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

// A Record contains alignment data for one BAM alignment record.
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

// RefID returns the target ID number for the alignment.
func (self *Record) RefID() int {
	self.unmarshalData()
	return int(self.tid())
}

// Name returns the name of the alignment query.
func (self *Record) Name() string {
	self.unmarshalData()
	return self.nameStr
}

// Seq returns a byte slice containing the sequence of the alignment query.
func (self *Record) Seq() []byte {
	self.unmarshalData()
	return self.seqBytes
}

// Quality returns an int8 slice containing the Phred quality scores of the alignment query.
func (self *Record) Quality() []int8 {
	self.unmarshalData()
	return self.qualScores
}

// Cigar returns a slice of CigarOps describing the alignment.
func (self *Record) Cigar() []CigarOp {
	self.unmarshalData()
	return self.cigar
}

// Tag returns an Aux tag whose tag ID matches the first two bytes of tag and true.
// If no tag matches, nil and false are returned.
func (self *Record) Tag(tag []byte) (v Aux, ok bool) {
	self.unmarshalData()
	for i := range self.auxTags {
		if bytes.Compare(self.auxTags[i][:2], tag) == 0 {
			return self.auxTags[i], true
		}
	}
	return
}

// Tags returns all Aux tags for the aligment.
func (self *Record) Tags() []Aux {
	self.unmarshalData()
	return self.auxTags
}

// Start returns the lower-coordinate end of the alignment.
func (self *Record) Start() int {
	return int(self.pos())
}

// Len returns the length of the alignment.
func (self *Record) Len() int {
	return int(self.lQseq())
}

// End returns the higher-coordinate end of the alignment.
// This is the start plus the sum of CigarMatch lengths.
func (self *Record) End() int {
	var mlen int
	for _, co := range self.Cigar() {
		if co.Type() == CigarMatch {
			mlen += co.Len()
		}
	}
	return int(self.pos()) + mlen
}

// Score returns the quality of the alignment.
func (self *Record) Score() byte {
	return self.qual()
}

// Flags returns the SAM flags for the alignment record.
func (self *Record) Flags() Flags {
	return self.flag()
}

// Strand returns an int8 indicating the strand of the alignment. A positive return indicates
// alignment in the forward orientation, a negative returns indicates alignemnt in the reverse
// orientation.
func (self *Record) Strand() int8 {
	if self.Flags()&Reverse == Reverse {
		return -1
	}
	return 1
}

// NextRefID returns the reference ID of the next segment/mate.
func (self *Record) NextRefID() int {
	return int(self.mtid())
}

// NextStart returns the start position of the next segment/mate.
func (self *Record) NextStart() int {
	return int(self.mpos())
}

// String returns a string representation of the Record.
func (self *Record) String() string {
	return fmt.Sprintf("%s %v %d:%d..%d %d %v %d:%d %d %s %v %v",
		self.Name(),
		self.Flags(),
		self.RefID(),
		self.Start(),
		self.End(),
		self.Score(),
		self.Cigar(),
		self.NextRefID(),
		self.NextStart(),
		self.Len(),
		self.Seq(),
		self.Quality(),
		self.Tags())
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

// unmarshalData interogates the bam1_t->data in the context of the bam1_t description fields to fill the Record's fields.
// unmarshalData is idempotent in this implementation although this may change.
func (self *Record) unmarshalData() {
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

// A CigarOp represents a Compact Idiosyncratic Gapped Alignment Report operation.
type CigarOp uint32

// Type returns the type of the CIGAR operation for the CigarOp.
func (co CigarOp) Type() CigarOpType { return CigarOpType(co & 0xf) }

// Len returns the number of positions affected by the CigarOp CIGAR operation.
func (co CigarOp) Len() int { return int(co >> 4) }

// String returns the string representation of the CigarOp
func (co CigarOp) String() string { return fmt.Sprintf("%d%s", co.Len(), co.Type().String()) }

// A CigarOpType represents the type of operation described by a CigarOp.
type CigarOpType byte

const (
	CigarMatch       CigarOpType = iota // Alignment match (can be a sequence match or mismatch).
	CigarInsertion                      // Insertion to the reference.
	CigarDeletion                       // Deletion from the reference.
	CigarSkipped                        // Skipped region from the reference.
	CigarSoftClipped                    // Soft clipping (clipped sequences present in SEQ).
	CigarHardClipped                    // Hard clipping (clipped sequences NOT present in SEQ).
	CigarPadded                         // Padding (silent deletion from padded reference).
	CigarEqual                          // Sequence match.
	CigarMismatch                       // Sequence mismatch.
	lastCigar
)

var cigarOps = []string{"M", "I", "D", "N", "S", "H", "P", "=", "X", "?"}

// String returns the string representation of a CigarOpType.
func (ct CigarOpType) String() string {
	if ct < 0 || ct > lastCigar {
		ct = lastCigar
	}
	return cigarOps[ct]
}

// An Aux represents an auxilliary tag data field from a SAM alignment record.
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

// parseAux examines the data of a SAM record's OPT fields,
// returning a slice of Aux that are backed by the original data.
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

// buildAux constructs a single byte slice that represents a slice of Aux.
func buildAux(aa []Aux) (aux []byte) {
	for _, a := range aa {
		// TODO: validate each 'a'
		// TODO: note that Z and H types should have a terminal null added.
		aux = append(aux, []byte(a)...)
	}
	return
}

// String returns the string representation of an Aux type.
func (self Aux) String() string {
	return fmt.Sprintf("%s:%c:%v", []byte(self[:2]), auxTypes[self.Type()], self.Value())
}

// A Tag represents an auxilliary tag label.
type Tag [2]byte

// String returns a string representation of a Tag.
func (self Tag) String() string { return string(self[:]) }

// Tag returns the string representation of the tag ID.
func (self Aux) Tag() Tag { var t Tag; copy(t[:], self[:2]); return t }

// Type returns a byte corresponding to the type of the auxilliary tag.
// Returned values are in {'A', 'c', 'C', 's', 'S', 'i', 'I', 'f', 'Z', 'H', 'B'}.
func (self Aux) Type() byte { return self[2] }

// Value returns v containing the value of the auxilliary tag.
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
