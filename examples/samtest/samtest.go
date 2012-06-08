package main

import (
	"fmt"
	"github.com/kortschak/boom"
)

func main() {
	sf, err := boom.OpenSAM("test.sam", "")
	if err != nil {
		panic(err)
	}
	for {
		r, _, err := sf.Read()
		if err != nil {
			break
		}
		fmt.Printf("%s %d-%d (%d) %v %d %d %s %v Fl:%#0.2x %q\n", r.Name(), r.Start(), r.End(), r.Len(), r.Cigar(), r.Strand(), r.Score(), r.Seq(), r.Quality(), r.Flags(), r.Tags())
	}
}
