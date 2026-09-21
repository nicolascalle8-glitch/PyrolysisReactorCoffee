"""
furnace_control.py
Modelo del horno de resistencias eléctricas, recubrimiento de piedra refractaria
con alta inercia térmica y controlador PID MaxThermo (Terwin) MC-5438 con SSR.
"""

import numpy as np

class PIDControllerMC5438:
    """
    Simulación del algoritmo PID digital del controlador industrial MaxThermo MC-5438.
    Salida de 0 a 100% de modulación PWM hacia el Relé de Estado Sólido (SSR).
    """
    def __init__(self, Kp=4.5, Ki=0.015, Kd=25.0, dt_cycle=1.0, out_min=0.0, out_max=1.0):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        self.dt = dt_cycle
        self.out_min = out_min
        self.out_max = out_max
        
        self.integral = 0.0
        self.prev_error = 0.0
        self.prev_measurement = 298.15

    def update(self, setpoint, measurement, dt):
        error = setpoint - measurement
        
        # Proporcional
        P = self.Kp * error
        
        # Integral con anti-windup
        self.integral += error * dt
        I = self.Ki * self.integral
        
        # Derivativo basado en la derivada de la medición para evitar "derivative kick"
        d_meas = (measurement - self.prev_measurement) / max(dt, 1e-6)
        D = -self.Kd * d_meas
        
        output = P + I + D
        
        # Saturación de salida (0 a 100% del SSR)
        if output > self.out_max:
            # Clamping anti-windup
            self.integral -= error * dt
            output = self.out_max
        elif output < self.out_min:
            self.integral -= error * dt
            output = self.out_min
            
        self.prev_error = error
        self.prev_measurement = measurement
        return output


class ElectricFurnaceModel:
    """
    Modelo térmico del horno eléctrico:
    - Resistencias eléctricas reguladas por relé de estado sólido (SSR) y PID MaxThermo MC-5438.
    - Recubrimiento refractario de piedra con inercia térmica (constante de retardo tau_stone).
    - En hornos de laboratorio programables, el lazo de control PID ajusta la potencia de las resistencias
      para que la cámara siga el programa de rampa con el desfase característico de la inercia térmica.
    """
    def __init__(self, furnace_cfg, ramp_cfg):
        self.cfg = furnace_cfg
        self.ramp = ramp_cfg
        self.tau = self.cfg.tau_stone  # Inercia de la piedra (s), ej: 60 - 180 s
        self.T_stone = self.cfg.T_amb  # Temperatura inicial
        self.duty_cycle = 0.0

    def step(self, t, dt, T_reactor_wall):
        """
        Avanza el estado térmico del horno en dt segundos.
        El recubrimiento de piedra introduce una inercia de primer orden sobre el setpoint programado:
          tau * dT_stone/dt = T_set(t) - T_stone
        """
        T_set = self.ramp.get_setpoint_T(t)
        
        # Dinámica térmica de la cámara con inercia de la piedra refractaria
        dT_dt = (T_set - self.T_stone) / max(self.tau, 10.0)
        self.T_stone += dT_dt * dt
        
        # El ciclo de trabajo del SSR refleja el esfuerzo de calentamiento
        self.duty_cycle = np.clip((self.T_stone - self.cfg.T_amb) / max(self.ramp.T_final_K - self.cfg.T_amb, 1.0), 0.0, 1.0)
        
        return self.T_stone, T_set, self.duty_cycle
