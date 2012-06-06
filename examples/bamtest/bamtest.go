package main

import (
	"fmt"
	"github.com/kortschak/boom"
)

func main() {
	bf, err := boom.OpenBAM("test.bam")
	if err != nil {
		panic(err)
	}
	for {
		r, n, err := bf.Read()
		if err != nil {
			break
		}
		fmt.Println(n, r.Start(), r.End(), r.Len(), r.Quality(), r.Flags())
	}
}
