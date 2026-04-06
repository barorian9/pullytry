import serial
import time

PORT = "COM4"          # change if needed
BAUD = 115200
FILENAME = "31.csv"

ser = serial.Serial(PORT, BAUD, timeout=1)
time.sleep(2)  # let Arduino reset fully

with open(FILENAME, "w") as f:
    print("Waiting for Arduino...")

    started = False

    while True:
        line = ser.readline().decode(errors="ignore").strip()

        if not line:
            continue

        # Wait for sync signal
        if line == "#START":
            print("Started logging")
            started = True
            continue

        # Write header once
        if started and line.startswith("t_s"):
            f.write(line + "\n")
            f.flush()
            continue

        # Write numeric data
        if started and line[0].isdigit():
            f.write(line + "\n")
            f.flush()

ser.close()
