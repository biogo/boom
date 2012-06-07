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
		fmt.Printf("%s %d %d %d %d %d %x, %q\n", r.ID(), r.Start(), r.End(), r.Len(), r.Strand(), r.Quality(), r.Flags(), r.Aux())
	}
}
