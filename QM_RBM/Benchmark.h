#pragma once
#include <iostream>
#include <chrono>

class Timer {
public:
	Timer(std::string timerName, System* system) {
		m_startTime = std::chrono::high_resolution_clock::now();
		m_system = system;
		name = timerName;
	}
	~Timer() {
		Stop();
	}

	void Stop() {
		auto endTime = std::chrono::high_resolution_clock::now();
		auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_startTime).time_since_epoch().count();
		auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endTime).time_since_epoch().count();
		
		double elapsed = (end - start)*0.001;

		std::cout << "Elapsed time of " << name << " was "<<elapsed << " milliseconds \n";
		m_system->setElapsedTime(elapsed);
	}

private:
	std::chrono::time_point< std::chrono::high_resolution_clock> m_startTime;
	std::string name = "";
	System* m_system = nullptr;
};