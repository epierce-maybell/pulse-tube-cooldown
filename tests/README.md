# Rules for the test suite

1. Acceptance tests, located in `tests/acceptance/`, are written by humans and never edited by agents. Interface acceptance tests ensure that the public interfaces are sane, intuitive, and complete. End-to-end acceptance tests ensure that the entire system behaves physically over a full simulation. Component acceptance tests ensure that individual components follow known solutions, preserve invariants, or otherwise do what they should.
2. Unit tests should enforce behavior, not implementation. 
3. The agent that writes tests for a unit of code should not be the agent that writes the code. While writing library code, agents should not modify tests.
4. The test suite should remain economical. If it takes more than a few minutes to run, that is a sign that the tests have been designed poorly.
