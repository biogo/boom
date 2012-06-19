package main

import (
	"code.google.com/p/biogo.boom"
	"fmt"
)

var (
	chr   = "gi|166362741|ref|NC_010296.1|"
	start = 3598704
	end   = 3598705
)

func main() {
	boom.Verbosity(0)
	bf, err := boom.OpenBAM("test-sort.bam")
	if err != nil {
		panic(err)
	}
	fmt.Printf("%s %d\n%s", bf.ReferenceNames(), bf.ReferenceLengths(), bf.Text())
	for {
		r, _, err := bf.Read()
		if err != nil {
			break
		}
		fmt.Printf("%d %s %d-%d (%d) %v %d %d %s %v Fl:%v %q\n", r.ReferenceID(), r.Name(), r.Start(), r.End(), r.Len(), r.Cigar(), r.Strand(), r.Score(), r.Seq(), r.Quality(), r.Flags(), r.Tags())
	}

	fmt.Println("Fetch:", chr, start, end)

	fn := func(r *boom.Record) {
		fmt.Printf("%d %s %d-%d (%d) %v %d %d %s %v Fl:%v %q\n", r.ReferenceID(), r.Name(), r.Start(), r.End(), r.Len(), r.Cigar(), r.Strand(), r.Score(), r.Seq(), r.Quality(), r.Flags(), r.Tags())
	}

	i, err := boom.LoadIndex("test-sort.bam")
	if err != nil {
		panic(err)
	}

	if tid, ok := bf.GetTargetID(chr); ok {
		_, err = bf.Fetch(i, tid, start, end, fn)
		if err != nil {
			panic(err)
		}
	} else {
		fmt.Printf("%q not present\n", chr)
	}
}
