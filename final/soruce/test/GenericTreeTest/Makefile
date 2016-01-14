generator: main.cpp
	g++ -O2 -std=c++11 main.cpp -DRECTEST -o rec.out
	g++ -O2 -std=c++11 main.cpp -DITETEST -o loop.out
testdata:
	g++ pin.cpp -std=c++11 -o pin.out
	./pin.out
test:
	time ./rec.out <in.txt
	time ./loop.out <in.txt
clean:
	rm *.out
	rm in.txt