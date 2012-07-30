package main

import (
	"code.google.com/p/biogo.boom"
	"fmt"
	"os"
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
	so, err := boom.OpenSAMFile(os.Stderr, "wh", bf.Header())
	if err != nil {
		panic(err)
	}
	defer so.Close()
	fmt.Printf("%s %d\n%s", bf.RefNames(), bf.RefLengths(), bf.Text())
	for {
		r, _, err := bf.Read()
		if err != nil {
			break
		}

		s := r.Seq()
		for i := range s {
			s[i] = 'S'
		}
		r.SetSeq(s)

		fmt.Printf("%d %s %d-%d (%d) %v %d %d %s %v Fl:%v %q\n", r.RefID(), r.Name(), r.Start(), r.End(), r.Len(), r.Cigar(), r.Strand(), r.Score(), r.Seq(), r.Quality(), r.Flags(), r.Tags())
		_, err = so.Write(r)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
		}
		fmt.Println()
	}

	fmt.Println("Fetch:", chr, start, end)

	fn := func(r *boom.Record) {
		fmt.Printf("%d %s %d-%d (%d) %v %d %d %s %v Fl:%v %q\n", r.RefID(), r.Name(), r.Start(), r.End(), r.Len(), r.Cigar(), r.Strand(), r.Score(), r.Seq(), r.Quality(), r.Flags(), r.Tags())
	}

	i, err := boom.LoadIndex("test-sort.bam")
	if err != nil {
		panic(err)
	}

	if tid, ok := bf.RefID(chr); ok {
		_, err = bf.Fetch(i, tid, start, end, fn)
		if err != nil {
			panic(err)
		}
	} else {
		fmt.Printf("%q not present\n", chr)
	}
}
