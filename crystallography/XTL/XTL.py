import serial
import time
import random
import logging
import time
from rich import print
import plotly.graph_objects as go
import pandas as pd
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import datetime
import threading

class Graph:
    def __init__(self, filename):
        self.filename = filename
        self.fig = go.Figure()
        self.x_data = []
        self.y1_data = []
        self.y2_data = []
        self.server_thread = None
    def update(self):
        try:
            with open(self.filename + '.txt', 'r') as file:
                lines = file.readlines()
                self.x_data = []
                self.y1_data = []
                self.y2_data = []
                for line in lines:
                    parts = line.strip().split(' - ')
                    if len(parts) == 2:
                        time_str, temps_str = parts
                        time_obj = datetime.datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S,%f')
                        time = (time_obj.hour * 3600) + (time_obj.minute * 60) + time_obj.second + (time_obj.microsecond / 1000000)
                        temps = temps_str.split(', ')
                        if len(temps) == 2:
                            temp1 = float(temps[0])
                            temp2 = float(temps[1])
                            self.x_data.append(time)
                            self.y1_data.append(temp1)
                            self.y2_data.append(temp2)
        except FileNotFoundError:
            pass
        self.fig.data = []
        self.fig.add_trace(go.Scatter(x=self.x_data, y=self.y1_data, name='Actual Temp'))
        self.fig.add_trace(go.Scatter(x=self.x_data, y=self.y2_data, name='Set Temp'))
    def show(self):
        def run_app():
            app = dash.Dash(__name__)
            app.layout = html.Div([
                dcc.Graph(id='live-graph', figure=self.fig),
                dcc.Interval(
                    id='graph-update',
                    interval=10000, # in milliseconds
                    n_intervals=0
                )
            ])
            @app.callback(
                Output('live-graph', 'figure'),
                [Input('graph-update', 'n_intervals')]
            )
            def update_graph_scatter(n):
                self.update()
                return self.fig
            app.run_server()
        self.server_thread = threading.Thread(target=run_app)
        self.server_thread.daemon = True  # Set as daemon so it stops when main program stops
        self.server_thread.start()
    def stop(self):
        if self.server_thread is not None:
            try:
                import psutil
                import os
                current_process = psutil.Process(os.getpid())
                for child in current_process.children(recursive=True):
                    child.terminate()
            except Exception as e:
                print(f"Error stopping graphing thread: {e}")

class MotorController:
    def __init__(self, port):
        self.port = port
        self.serial_connection = serial.Serial(port, 115200, timeout=2)
        self.motor_addresses = self.get_motor_addresses()
        self.motors = [Motor(self.serial_connection, address) for address in self.motor_addresses]
    def get_motor_addresses(self):
        motor_addresses = []
        for i in range(1, 100):  # Assume max 99 devices
            command = f'/{i} get comm.address\r\n'
            self.serial_connection.write(command.encode())
            response = self.serial_connection.readline().decode()
            if '@' in response:
                motor_addresses.append(i)
            else:
                break
        return motor_addresses
    def create_motor(self, address):
        return Motor(self.serial_connection, address)
    def set_rotation_speed(self, speed):
        for motor in self.motors:
            motor.set_speed(speed)

# Module 1: Motor Control
class Motor:
    def __init__(self, serial_connection, address):
        self.serial_connection = serial_connection
        self.address = address
        self.speed = 0
        self.direction_multiplier = random.choice([-1, 1])
        self.last_change_time = time.time()
        self.next_change_duration = random.uniform(5, 10)
    def set_speed(self, speed):
        self.speed = speed  # Add this line to update the speed attribute
        command = f"/{self.address} set maxspeed {speed}\r\n"
        # print(f"Sending command to motor {self.address}: {command}")
        self.serial_connection.write(command.encode())
        time.sleep(0.1)  # Add a 100ms delay
        # response = self.serial_connection.readline()
        # print(f"Received response from motor {self.address}: {response}")
        # if b'OK' in response:
        #     print("Response OK")
        # else:
        #     print("Error: Invalid response received from motor controller")
    def move_relative(self, distance):
        command = f"/{self.address} move rel {distance}\r\n"
        self.serial_connection.write(command.encode())
        time.sleep(0.1)  # Add a 100ms delay
        # response = self.serial_connection.readline()
        # if response.startswith(f'@{self.address:02}'.encode()):
        #     print("Response OK")
        # else:
        #     print("Error: Invalid response received from motor controller")
    def update_direction(self):
        current_time = time.time()
        if current_time - self.last_change_time >= self.next_change_duration:
            self.direction_multiplier *= -1
            self.last_change_time = current_time
            self.next_change_duration = random.uniform(5, 10)
            # print(f"Motor {self.address} direction changed. Next direction change will occur in {self.next_change_duration} seconds")
    def move(self):
        distance = int(self.speed * self.direction_multiplier)
        # print(f"Motor {self.address} speed: {self.speed}")
        # print(f"Motor {self.address} direction multiplier: {self.direction_multiplier}")
        # print(f"Motor {self.address} distance: {distance}")
        if distance == 0:
            print(f"Motor {self.address} has zero distance to move. Check speed and direction.")
        command = f"/{self.address} move rel {distance}\r\n"
        # print(f"Sending command to motor {self.address}: {command}")
        self.serial_connection.write(command.encode())
        time.sleep(0.1)  # Add a 100ms delay
        response = self.serial_connection.readline()
        # print(f"Received response from motor {self.address}: {response}")
        # if b'OK' in response:
        #     print("Response OK")

# Module 2: Heater Control
class Heater:
    def __init__(self, port):
        self.port = port
        self.serial_connection = serial.Serial(port, 9600, timeout=1)

    def read_temperature(self):
        command = b'IN_PV_2'
        self.serial_connection.write(command + b'\n')
        time.sleep(0.1)  # Add a 100ms delay
        response = self.serial_connection.readline().decode().strip()
        return response

    def set_temperature(self, temperature):
        command = f"OUT_SP_1 {temperature}\r\n"
        self.serial_connection.write(command.encode())
        time.sleep(0.1)  # Add a 100ms delay
        # print(f"Temperature set {temperature}")

    def start_tempering(self):
        command = b'START_1'
        self.serial_connection.write(command + b'\n')
        time.sleep(0.1)  # Add a 100ms delay
        print("Tempering function started")

    def stop_tempering(self):
        command = b'STOP_1'
        self.serial_connection.write(command + b'\n')
        time.sleep(0.1)  # Add a 100ms delay
        print("Tempering function stopped")

# Module 3: Temperature Profile
class TemperatureProfile:
    def __init__(self, filename):
        self.filename = filename
        self.profile = []
        self.start_time = None
    def load_profile(self):
        with open(self.filename, 'r') as file:
            for line in file:
                line = line.strip().replace('\n', '')  # Remove newline characters
                parts = line.split(' -- ')
                if len(parts) == 2:  # Check if the line was split into two parts
                    try:
                        temperature, time = float(parts[0]), float(parts[1])
                        self.profile.append((temperature, time))
                    except ValueError:
                        print(f"Invalid value in line: {line}")
    def get_current_temperature(self, elapsed_time):
        total_time = sum(time for _, time in self.profile)
        if elapsed_time >= total_time * 3600:
            return self.profile[-1][0]
        
        current_step = 0
        cumulative_time = 0
        previous_temp = 26 # assuming starting temperature is 26 C
        for i, (temperature, time) in enumerate(self.profile):
            step_time = time * 3600 # convert hours to seconds
            if elapsed_time <= cumulative_time + step_time:
                fraction = (elapsed_time - cumulative_time) / step_time
                current_temp = previous_temp + (temperature - previous_temp) * fraction
                return current_temp
            previous_temp = temperature
            cumulative_time += step_time
    def update_temperature(self, time):
        current_temperature = self.get_current_temperature(time)
        return current_temperature

# Module 4: Main Program
class Program:
    def __init__(self):
        self.port = "COM5"
        self.controller = MotorController(self.port)
        self.heater = Heater("COM4")
        self.profile = TemperatureProfile("profile.txt")
        self.filename = 'temperature_log'
        self.graph_initialized = False  # Add this line

    def display_menu(self):
        print("Menu:")
        print("1. Set file name for temp data:")
        print("2. Set rotation speed")
        print("3. Start and hold heater")
        print("4. Start program")
        print("5. Cancel program and exit")

    def handle_input(self, choice):
        if choice == "3":
            self.start_and_hold_heater()
        elif choice == "4":
            self.initialize_graph()  # Start the graph here
            if all(motor.speed != 0 for motor in self.controller.motors):
                self.start_program()
            else:
                print("Speed not set")
                self.display_menu()
                input_choice = input("Enter choice: ")
                self.handle_input(input_choice)
        elif choice == "5":
            self.cancel_program_and_exit()
        elif choice == "2":
            speed = int(input("Enter rotation speed: "))
            self.controller.set_rotation_speed(speed*10000)
        elif choice == "1":
            filename = input("Enter file name for temp data: ")
            global logger
            logger.handlers.clear()  # Remove existing handlers
            handler = TimedFileHandler(filename + '.txt')
            handler.setLevel(logging.INFO)
            formatter = logging.Formatter('%(asctime)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            self.filename = filename
            # self.initialize_graph()

    def start_and_hold_heater(self):
        temperature = float(input("Enter temperature: "))
        self.heater.set_temperature(temperature)
        self.heater.start_tempering()

    def start_program(self):
        self.profile.load_profile()
        start_time = time.time()
        while True:
            current_time = time.time()
            elapsed_time = current_time - start_time
            current_temperature = self.profile.get_current_temperature(elapsed_time)
            actual_temp = float(self.heater.read_temperature().split(' ')[0])
            logger.info(f'{actual_temp}, {current_temperature}')
            # print(f"Current temperature: {current_temperature}")
            for motor in self.controller.motors:
                motor.update_direction()
                motor.move()
            self.heater.set_temperature(current_temperature)
            time.sleep(1)

    def cancel_program_and_exit(self):
        self.heater.stop_tempering()
        sys.exit(0)

    def set_rotation_speed(self):
        speed = int(input("Enter rotation speed (microsteps per second): "))
        for motor in self.motors:
            motor.set_speed(speed)

    def cancel_program_and_exit(self):
        if self.graph_initialized:
            self.graph.stop()  # Stop the graphing thread
        self.heater.stop_tempering()
        sys.exit(0)
    def initialize_graph(self):
        self.graph = Graph(self.filename)  # Store the graph instance
        self.graph.show()
        self.graph_initialized = True

class TimedFileHandler(logging.FileHandler):
    def __init__(self, filename, mode='a', encoding=None, delay=False):
        super().__init__(filename, mode, encoding, delay)
        self.last_write_time = time.time()
    def emit(self, record):
        current_time = time.time()
        if current_time - self.last_write_time >= 10:
            self.last_write_time = current_time
            super().emit(record)

# Create a logger
logger = logging.getLogger('temperature_logger')
logger.setLevel(logging.INFO)
# Create a file handler
handler = TimedFileHandler('temperature_log.txt')
handler.setLevel(logging.INFO)
# Create a formatter and set it for the handler
formatter = logging.Formatter('%(asctime)s - %(message)s')
handler.setFormatter(formatter)
# Add the handler to the logger
logger.addHandler(handler)

def main():
    program = Program()
    while True:
        program.display_menu()
        input_choice = input("Enter choice: ")
        program.handle_input(input_choice)

if __name__ == "__main__":
    main()