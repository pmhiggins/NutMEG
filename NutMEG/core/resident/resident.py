from abc import ABC, abstractmethod


class Resident(ABC):
    """
    Abstract class for an ecosystem property that evolves in time.
    This enforces that any Resident has the ability to take a numerical
    time-step.
    """

    @abstractmethod
    def take_step(self, dt):
        """ Advance resident by time-step dt. """
        pass
